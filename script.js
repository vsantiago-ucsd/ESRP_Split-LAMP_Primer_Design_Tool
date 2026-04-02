// Global state object - single source of truth
let designState = {
    architecture: 'f2-only',
    inputs: {
        mirna1: { name: '', sequence: '', parsed: '', length: 0, gc: 0 },
        mirna2: { name: '', sequence: '', parsed: '', length: 0, gc: 0 }
    },
    options: {
        generateLoop: true,
        optimizeSpecificity: true,
        checkHairpins: true,
        calculateDG: true
    },
    outputs: {
        template: { seq: '', len: 0, gc: 0 },
        lf: { seq: '', len: 0, gc: 0, tm: 0, dg: 0, hairpin: 'None' },
        lb: { seq: '', len: 0, gc: 0, tm: 0, dg: 0, hairpin: 'None' },
        fip: { seq: '', len: 0, gc: 0, tm: 0, dg: 0, hairpin: 'None' },
        bip: { seq: '', len: 0, gc: 0, tm: 0, dg: 0, hairpin: 'None' },
        f2: { seq: '', len: 0, gc: 0, tm: 0, dg: 0 },
        b2: { seq: '', len: 0, gc: 0, tm: 0, dg: 0 },
        f1c: { seq: '', len: 0, gc: 0, tm: 0, dg: 0 },
        b1c: { seq: '', len: 0, gc: 0, tm: 0, dg: 0 }
    },
    // Snapshot written once at generation — never mutated, used by Reset
    original: {
        f2: '',
        b2: '', 
        f1c: '', 
        b1c: '',
        lf: '',
        lb: ''
    }
};

const NN_PARAMS = {
    'AA': { dH: -7.9, dS: -22.2 }, 'TT': { dH: -7.9, dS: -22.2 },
    'AT': { dH: -7.2, dS: -20.4 }, 'TA': { dH: -7.2, dS: -21.3 },
    'CA': { dH: -8.5, dS: -22.7 }, 'TG': { dH: -8.5, dS: -22.7 },
    'GT': { dH: -8.4, dS: -22.4 }, 'AC': { dH: -8.4, dS: -22.4 },
    'CT': { dH: -7.8, dS: -21.0 }, 'AG': { dH: -7.8, dS: -21.0 },
    'GA': { dH: -8.2, dS: -22.2 }, 'TC': { dH: -8.2, dS: -22.2 },
    'CG': { dH: -10.6, dS: -27.2 }, 'GC': { dH: -9.8, dS: -24.4 },
    'GG': { dH: -8.0, dS: -19.9 }, 'CC': { dH: -8.0, dS: -19.9 }
};

// Normalize sequence parsing - single function used everywhere
function parseSequence(input) {
    // Remove FASTA header(s)
    let cleaned = input.replace(/^>.*$/gm, '');
    
    // Strip whitespace, numbers, and line breaks
    cleaned = cleaned.replace(/[^a-zA-Z]/g, '');
    
    // Convert to uppercase
    cleaned = cleaned.toUpperCase();
    
    // Validate alphabet - ACGTU and T (DNA) all allowed
    const validBases = /^[ACGTU]+$/;
    if (cleaned && !validBases.test(cleaned)) {
        return { valid: false, sequence: '', error: 'Invalid characters. Only A, C, G, T/U allowed.' };
    }
    
    // Convert U to T for DNA oligos
    const dnaSequence = cleaned.replace(/U/g, 'T');
    
    return { valid: true, sequence: dnaSequence, error: '' };
}

// Calculate GC content
function calculateGC(sequence) {
    if (!sequence) return 0;
    const gcCount = (sequence.match(/[GC]/g) || []).length;
    return Math.round((gcCount / sequence.length) * 100);
}

// Strictly SantaLucia 1998 compliant thermodynamic calculation
function calcThermo(sequence) {
    sequence = sequence.toUpperCase();

    // 1. Sum nearest-neighbor parameters
    let totalDH = 0;
    let totalDS = 0;

    for (let i = 0; i < sequence.length - 1; i++) {
        const dinuc = sequence.substring(i, i + 2);
        if (NN_PARAMS[dinuc]) {
            totalDH += NN_PARAMS[dinuc].dH;
            totalDS += NN_PARAMS[dinuc].dS;
        }
    }

    // 2. Initiation — one correction per terminal end (SantaLucia 1998, Table 2)
    const fiveEnd  = sequence[0];
    const threeEnd = sequence[sequence.length - 1];

    // 5' terminal
    if (fiveEnd === 'G' || fiveEnd === 'C') {
        totalDH += 0.1;
        totalDS += -2.8;
    } else {
        totalDH += 2.3;
        totalDS += 4.1;
    }

    // 3' terminal
    if (threeEnd === 'G' || threeEnd === 'C') {
        totalDH += 0.1;
        totalDS += -2.8;
    } else {
        totalDH += 2.3;
        totalDS += 4.1;
    }

    // 3. Symmetry correction (only if self-complementary)
    if (sequence === reverseComplement(sequence)) {
        totalDS -= 1.4;
    }

    return { dH: totalDH, dS: totalDS };
}

function calculateTm(sequence, naConc = 50, mgConc = 8, primerConc = 0.25, dntpConc = 0.8) {
    if (!sequence) return 0;
    sequence = sequence.toUpperCase();

    // Short sequence fallback (Wallace rule) — SantaLucia recommends NN for all lengths
    // but <8 bp is unreliable with NN; keeping your threshold at 14 is reasonable
    if (sequence.length < 14) {
        const a = (sequence.match(/A/g) || []).length;
        const t = (sequence.match(/T/g) || []).length;
        const g = (sequence.match(/G/g) || []).length;
        const c = (sequence.match(/C/g) || []).length;
        return 2 * (a + t) + 4 * (g + c);
    }

    const { dH, dS } = calcThermo(sequence);

    // Salt correction on ΔS (SantaLucia 1998 Eq. 6 / Owczarzy approximation)
    const effectiveMg  = Math.max(0, mgConc - dntpConc);
    const naEffective  = naConc + 120 * Math.sqrt(effectiveMg);  // mM equivalent
    const saltCorr     = 0.368 * (sequence.length - 1) * Math.log(naEffective / 1000);
    const correctedDS  = dS + saltCorr;

    // Tm equation (SantaLucia 1998, Eq. 3) — non-self-complementary: CT/4
    const R  = 1.987;       // cal/mol·K
    const Ct = primerConc / 1e6;  // convert µM → M
    const Tm = (dH * 1000) / (correctedDS + R * Math.log(Ct / 4)) - 273.15;

    return Math.round(Tm * 10) / 10;
}

function calculateDeltaG(sequence, temperature = 37, naConc = 50, mgConc = 8, dntpConc = 0.8) {
    if (!sequence) return 0;
    sequence = sequence.toUpperCase();

    const { dH, dS } = calcThermo(sequence);

    // Salt correction on ΔS
    const effectiveMg = Math.max(0, mgConc - dntpConc);
    const naEffective = naConc + 120 * Math.sqrt(effectiveMg);
    const saltCorr    = 0.368 * (sequence.length - 1) * Math.log(naEffective / 1000);
    const correctedDS = dS + saltCorr;

    // ΔG at reaction temperature (65°C for LAMP)
    const tempK  = temperature + 273.15;
    const deltaG = dH - (tempK * correctedDS / 1000);

    return Math.round(deltaG * 100) / 100;
}

// Calculate ΔG for the 5' end (first 6 bp)
function calculateDeltaG5Prime(sequence, temperature = 37, naConc = 50, mgConc = 8, dntpConc = 0.8) {
    if (!sequence) return 0;
    let fiveEnd = sequence.length > 6 ? sequence.slice(0,6) : sequence;
    return calculateDeltaG(fiveEnd, temperature, naConc, mgConc, dntpConc);
}

// Calculate ΔG for the 3' end (last 6 bp)
function calculateDeltaG3Prime(sequence, temperature = 37, naConc = 50, mgConc = 8, dntpConc = 0.8) {
    if (!sequence) return 0;
    let threeEnd = sequence.length > 6 ? sequence.slice(-6) : sequence;
    return calculateDeltaG(threeEnd, temperature, naConc, mgConc, dntpConc);
}

// Parse and validate sequence with immediate feedback
function parseAndValidateSequence(mirnaId) {
    const textarea = document.getElementById(`${mirnaId}-sequence`);
    const errorDiv = document.getElementById(`${mirnaId}-error`);
    const infoDiv = document.getElementById(`${mirnaId}-info`);
    const lengthSpan = document.getElementById(`${mirnaId}-length`);
    const gcSpan = document.getElementById(`${mirnaId}-gc`);
    
    const result = parseSequence(textarea.value);
    
    if (!textarea.value) {
        errorDiv.classList.remove('visible');
        infoDiv.classList.remove('visible');
        return;
    }
    
    if (!result.valid) {
        errorDiv.textContent = result.error;
        errorDiv.classList.add('visible');
        infoDiv.classList.remove('visible');
        return;
    }
    
    errorDiv.classList.remove('visible');
    
    // Calculate and show length and GC
    const length = result.sequence.length;
    const gc = calculateGC(result.sequence);
    
    lengthSpan.textContent = length;
    gcSpan.textContent = gc;
    infoDiv.classList.add('visible');
    
    // Update state
    if (mirnaId === 'mirna1') {
        designState.inputs.mirna1.parsed = result.sequence;
        designState.inputs.mirna1.length = length;
        designState.inputs.mirna1.gc = gc;
    } else {
        designState.inputs.mirna2.parsed = result.sequence;
        designState.inputs.mirna2.length = length;
        designState.inputs.mirna2.gc = gc;
    }
}

// Tab switching
document.querySelectorAll('.tab').forEach(tab => {
    tab.addEventListener('click', function() {
        const targetTab = this.dataset.tab;
        
        document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
        this.classList.add('active');
        
        document.querySelectorAll('.tab-content').forEach(content => content.classList.remove('active'));
        document.getElementById(targetTab).classList.add('active');
    });
});

// Architecture selection
function selectArchitecture(arch) {
    designState.architecture = arch;
    
    // Update radio selection
    document.querySelectorAll('input[name="architecture"]').forEach(radio => {
        radio.checked = (radio.value === arch);
    });
    
    // Update visual selection
    document.querySelectorAll('.architecture-option').forEach(option => {
        option.classList.remove('selected');
    });
    event.currentTarget.classList.add('selected');
    
    // Show/hide second microRNA section
    const mirna2Section = document.getElementById('mirna2-section');
    if (arch === 'f2-and-b2') {
        mirna2Section.classList.add('visible');
    } else {
        mirna2Section.classList.remove('visible');
    }
}

// Set state - idle, loading, or results
function setState(state) {
    const idleState = document.getElementById('idle-state');
    const loadingState = document.getElementById('loading-state');
    const resultsState = document.getElementById('results-state');
    
    idleState.classList.remove('hidden');
    loadingState.classList.remove('visible');
    resultsState.classList.remove('visible');
    
    if (state === 'loading') {
        idleState.classList.add('hidden');
        loadingState.classList.add('visible');
    } else if (state === 'results') {
        idleState.classList.add('hidden');
        resultsState.classList.add('visible');
    }
}

// Generate primers function
function generatePrimers() {
    // Validate inputs
    const mirna1Name = document.getElementById('mirna1-name').value.trim();
    const mirna1Seq = document.getElementById('mirna1-sequence').value.trim();
    
    if (!mirna1Name || !mirna1Seq) {
        alert('Please enter the first microRNA name and sequence');
        return;
    }
    
    const mirna1Result = parseSequence(mirna1Seq);
    if (!mirna1Result.valid) {
        alert('Invalid first microRNA sequence: ' + mirna1Result.error);
        return;
    }
    
    // If F2-and-B2 architecture, validate second microRNA
    if (designState.architecture === 'f2-and-b2') {
        const mirna2Name = document.getElementById('mirna2-name').value.trim();
        const mirna2Seq = document.getElementById('mirna2-sequence').value.trim();
        
        if (!mirna2Name || !mirna2Seq) {
            alert('Please enter the second microRNA name and sequence for F2-and-B2 architecture');
            return;
        }
        
        const mirna2Result = parseSequence(mirna2Seq);
        if (!mirna2Result.valid) {
            alert('Invalid second microRNA sequence: ' + mirna2Result.error);
            return;
        }
        
        // Update state
        designState.inputs.mirna2.name = mirna2Name;
        designState.inputs.mirna2.sequence = mirna2Seq;
        designState.inputs.mirna2.parsed = mirna2Result.sequence;
    }
    
    // Update state
    designState.inputs.mirna1.name = mirna1Name;
    designState.inputs.mirna1.sequence = mirna1Seq;
    designState.inputs.mirna1.parsed = mirna1Result.sequence;
    
    // Show loading
    setState('loading');
    document.querySelector('[data-tab="output"]').click();
    
    setTimeout(() => {
        
        // Update target info listed at top
        document.getElementById('result-architecture').textContent = 
            designState.architecture === 'f2-only' ? 'F2 Only' : 'F2 and B2 (AND-gate)';
        
        document.getElementById('result-target1-name').textContent = mirna1Name;
        const r1seqEl = document.getElementById('result-target1-seq');
        if (r1seqEl) r1seqEl.textContent = mirna1Result.sequence;
        
        if (designState.architecture === 'f2-and-b2') {
            const mirna2Result = parseSequence(document.getElementById('mirna2-sequence').value);
            document.getElementById('result-target2-info').style.display = 'block';
            document.getElementById('result-target2-name').textContent = 
                document.getElementById('mirna2-name').value;
            const r2seqEl = document.getElementById('result-target2-seq');
            if (r2seqEl) r2seqEl.textContent = mirna2Result.sequence;
        } else {
            document.getElementById('result-target2-info').style.display = 'none';
        }
        
        // Mock primer generation - replace with actual algorithm
        
        // Define all primer sequences
        const lfSeq = 'TCACTGATCTGGCCGTAGACCA';
        const lbSeq = 'TGACAGGACATCGGTGACAGT';
        const f1cSeq = 'CGGAGAGGTCGCGATAGTCA';
        const b1cSeq = 'GATGACAGTGACATCCTGCCT';

        let f2Seq = mirna1Result.sequence.substring(2, mirna1Result.sequence.length);
        let fipSeq = 'CGGAGAGGTCGCGATAGTCAT' + f2Seq; // f1c + T + f2
        let b2Seq;
        let bipSeq;

        // Save primers to designState for template cascade
        designState.outputs.lf.seq = lfSeq;
        designState.outputs.lb.seq = lbSeq;
        designState.outputs.f1c.seq = f1cSeq;
        designState.outputs.b1c.seq = b1cSeq;
        designState.outputs.f2.seq = f2Seq;
        designState.outputs.fip.seq = fipSeq;
        // Lock in the original — never overwritten, used by Reset
        designState.original.lf = lfSeq;
        designState.original.lb = lbSeq
        designState.original.f1c = f1cSeq;
        designState.original.b1c = b1cSeq;
        designState.original.f2 = f2Seq;
        
        // Loop primers (LF, LB)
        updatePrimerOutput('lf-seq', lfSeq, lfSeq.length, 
            calculateGC(lfSeq), 
            calculateTm(lfSeq), 
            calculateDeltaG5Prime(lfSeq), 
            calculateDeltaG3Prime(lfSeq)
            );
        
        updatePrimerOutput('lb-seq', lbSeq, lbSeq.length, 
            calculateGC(lbSeq), 
            calculateTm(lbSeq), 
            calculateDeltaG5Prime(lbSeq), 
            calculateDeltaG3Prime(lbSeq)
            );
        
        // Inner primers (FIP with F2, F1C, B1C)
        updatePrimerOutput('fip-seq', fipSeq, fipSeq.length, 
            calculateGC(fipSeq), 
            undefined, 
            undefined, 
            undefined,
            );

        updatePrimerOutput('f2-seq', f2Seq, f2Seq.length, 
            calculateGC(f2Seq), 
            calculateTm(f2Seq), 
            calculateDeltaG5Prime(f2Seq), 
            calculateDeltaG3Prime(f2Seq));
        
        updatePrimerOutput('f1c-seq', f1cSeq, f1cSeq.length, 
            calculateGC(f1cSeq), 
            calculateTm(f1cSeq), 
            calculateDeltaG5Prime(f1cSeq), 
            calculateDeltaG3Prime(f1cSeq));
        
        updatePrimerOutput('b1c-seq', b1cSeq, b1cSeq.length, 
            calculateGC(b1cSeq), 
            calculateTm(b1cSeq), 
            calculateDeltaG5Prime(b1cSeq), 
            calculateDeltaG3Prime(b1cSeq));

        // Architecture-specific primers (BIP, B2)
        if (designState.architecture === 'f2-and-b2') {
            const mirna2Result = parseSequence(document.getElementById('mirna2-sequence').value);
            b2Seq = mirna2Result.sequence;
            bipSeq = 'GATGACAGTGACATCCTGCCTA' + b2Seq.substring(1); // b1c + A + (b2 - first T)
            
            // Save primers to designState for template cascade
            designState.outputs.b2.seq = b2Seq;
            designState.outputs.bip.seq = bipSeq;
            // Lock in the original — never overwritten, used by Reset
            designState.original.b2 = b2Seq;

            // Generate template

            const { rawSequence: templateRaw }
                = generateTemplateUltramer(fipSeq, lfSeq, f1cSeq, b1cSeq, lbSeq, bipSeq);
            designState.outputs.template.seq = templateRaw;
            
            // BIP and B2 for two-input architecture
            updatePrimerOutput('bip-seq', bipSeq, bipSeq.length, 
                calculateGC(bipSeq), 
                undefined, 
                undefined, 
                undefined
                );
            
            updatePrimerOutput('b2-seq', b2Seq, b2Seq.length, 
                calculateGC(b2Seq), 
                calculateTm(b2Seq), 
                calculateDeltaG5Prime(b2Seq), 
                calculateDeltaG3Prime(b2Seq));

            updatePrimerOutput("template-seq", templateRaw, templateRaw.length, calculateGC(templateRaw));
        } else {
            // Single-input architecture
            b2Seq = 'TGGCAGTGTCTTAGCTGGTTGT';
            bipSeq = 'GATGACAGTGACATCCTGCCTA' + b2Seq.substring(1);

            // Save primers to designState for template cascade
            designState.outputs.b2.seq = b2Seq;
            designState.outputs.bip.seq = bipSeq;
            // Lock in the original — never overwritten, used by Reset
            designState.original.b2 = b2Seq;
            
            const { rawSequence: templateRaw }
                = generateTemplateUltramer(fipSeq, lfSeq, f1cSeq, b1cSeq, lbSeq, bipSeq);
            designState.outputs.template.seq = templateRaw;
            updatePrimerOutput("template-seq", templateRaw, templateRaw.length, calculateGC(templateRaw));
            
            updatePrimerOutput('bip-seq', bipSeq, bipSeq.length, 
                calculateGC(bipSeq), 
                undefined, 
                undefined, 
                undefined
                );
            
            updatePrimerOutput('b2-seq', b2Seq, b2Seq.length, 
                calculateGC(b2Seq), 
                calculateTm(b2Seq), 
                calculateDeltaG5Prime(b2Seq), 
                calculateDeltaG3Prime(b2Seq));
        }
        
        // Show results
        setState('results');

        // Enable F2 input and attach live-edit listener (attach once per generation)
        const f2Input = document.getElementById('f2-seq');
        f2Input.disabled = false;
        document.getElementById('f2-reset').disabled = false;
        f2Input.removeEventListener('input', primerInputHandler('f2'));
        f2Input.addEventListener('input', primerInputHandler('f2'));

        // Enable B2 input and attach live-edit listener (attach once per generation)
        const b2Input = document.getElementById('b2-seq');
        b2Input.disabled = false;
        document.getElementById('b2-reset').disabled = false;
        b2Input.removeEventListener('input', primerInputHandler('b2'));
        b2Input.addEventListener('input', primerInputHandler('b2'));

        // Enable F1C input and attach live-edit listener (attach once per generation)
        const f1cInput = document.getElementById('f1c-seq');
        f1cInput.disabled = false;
        document.getElementById('f1c-reset').disabled = false;
        f1cInput.removeEventListener('input', primerInputHandler('f1c'));
        f1cInput.addEventListener('input', primerInputHandler('f1c'));

        // Enable B1C input and attach live-edit listener (attach once per generation)
        const b1cInput = document.getElementById('b1c-seq');
        b1cInput.disabled = false;
        document.getElementById('b1c-reset').disabled = false;
        b1cInput.removeEventListener('input', primerInputHandler('b1c'));
        b1cInput.addEventListener('input', primerInputHandler('b1c'));

        // Enable LF input and attach live-edit listener (attach once per generation)
        const lfInput = document.getElementById('lf-seq');
        lfInput.disabled = false;
        document.getElementById('lf-reset').disabled = false;
        lfInput.removeEventListener('input', primerInputHandler('lf'));
        lfInput.addEventListener('input', primerInputHandler('lf'));

        // Enable LB input and attach live-edit listener (attach once per generation)
        const lbInput = document.getElementById('lb-seq');
        lbInput.disabled = false;
        document.getElementById('lb-reset').disabled = false;
        lbInput.removeEventListener('input', primerInputHandler('lb'));
        lbInput.addEventListener('input', primerInputHandler('lb'));
    updateDimerBox();}, 1500);
}

// Generic handler for sequence editing
function primerInputHandler(primerName) {
    return function() {
        onPrimerChange(primerName, this.value);
    };
}

// Generic live edit handler for any primer
function onPrimerChange(primerName, rawValue) {
    const errorDiv = document.getElementById(`${primerName}-edit-error`);
    const result = parseSequence(rawValue);

    if (!result.valid) {
        errorDiv.textContent = result.error;
        errorDiv.classList.add('visible');
        return;
    }
    errorDiv.classList.remove('visible');

    const sequence = result.sequence;

    // Update stats and designState
    designState.outputs[primerName].seq = sequence;
    updatePrimerOutput(
        `${primerName}-seq`, 
        sequence, 
        sequence.length, 
        calculateGC(sequence),
        calculateTm(sequence), 
        calculateDeltaG5Prime(sequence), 
        calculateDeltaG3Prime(sequence)
    );

    // Recompute dependent sequences based on primer type
    if (primerName === 'f2') {
        const f1cSeq = designState.outputs.f1c.seq;
        const newFip = f1cSeq + 'T' + sequence;
        designState.outputs.fip.seq = newFip;
        updatePrimerOutput('fip-seq', newFip, newFip.length, calculateGC(newFip),
            undefined, undefined, undefined);
    } else if (primerName === 'b2') {
        const b1cSeq = designState.outputs.b1c.seq;
        const newBip = b1cSeq + 'A' + sequence;
        designState.outputs.bip.seq = newBip;
        updatePrimerOutput('bip-seq', newBip, newBip.length, calculateGC(newBip),
            undefined, undefined, undefined);
    } else if (primerName == 'f1c'){
        const f2Seq = designState.outputs.f2.seq;
        const newFip = sequence + 'T' + f2Seq;
        designState.outputs.fip.seq = newFip;
        updatePrimerOutput('fip-seq', newFip, newFip.length, calculateGC(newFip),
            undefined, undefined, undefined);
    } else if (primerName == 'b1c'){
        const b2Seq = designState.outputs.b2.seq;
        const newBip = sequence + 'A' + b2Seq;
        designState.outputs.bip.seq = newBip;
        updatePrimerOutput('bip-seq', newBip, newBip.length, calculateGC(newBip),
            undefined, undefined, undefined);
    }

    // Recompute template ultramer
    const { rawSequence: newTemplate } = generateTemplateUltramer(
        designState.outputs.fip.seq,
        designState.outputs.lf.seq,
        designState.outputs.f1c.seq,
        designState.outputs.b1c.seq,
        designState.outputs.lb.seq,
        designState.outputs.bip.seq
    );
        designState.outputs.template.seq = newTemplate;

    // updatePrimerOutput handles innerHTML via highlightTemplate internally
    updatePrimerOutput('template-seq', newTemplate, newTemplate.length, calculateGC(newTemplate));

    updateDimerBox();
}

// Reset primer to originally generated sequence
function resetPrimer(primerName) {
    const original = designState.original[primerName];
    if (!original) return;

    const inputEl = document.getElementById(`${primerName}-seq`);
    inputEl.value = original;
    document.getElementById(`${primerName}-edit-error`).classList.remove('visible');

    // Re-run the normal change handler so dependent sequences (FIP, template, etc.) update
    onPrimerChange(primerName, original);
}

// Update primer output with stable IDs
function updatePrimerOutput(seqId, sequence, length, gc, tm, dg5, dg3) {
    const seqElement = document.getElementById(seqId);
    const primerPrefix = seqId.replace('-seq', '');
    
    // Update sequence — handle both <input> and <div> elements
        if (seqElement.tagName === 'INPUT') {
        seqElement.value = sequence;
    } else {
        if (primerPrefix === "template") {
            seqElement.innerHTML = highlightTemplate(sequence);
            // Update template hairpin stat here, after highlightTemplate runs
            const th = document.getElementById('template-hairpin');
            if (th) {
                const hr = detectHairpin(sequence);
                th.textContent = hr.summary;
                th.className = hr.summary === 'None' ? 'stat-value good'
                            : hr.strong.length > 0  ? 'stat-value bad'
                            : 'stat-value warning';
            }
        } else {
            seqElement.textContent = sequence;
        }
    }
    seqElement.setAttribute('data-sequence', sequence);
    
    // Enable copy button
    const copyBtn = seqElement.parentElement.querySelector('.btn-copy');
    if (copyBtn) copyBtn.disabled = false;
    
    // Update stats
    const lenEl = document.getElementById(`${primerPrefix}-len`);
    const gcEl  = document.getElementById(`${primerPrefix}-gc`);
    if (lenEl) lenEl.textContent = length + ' bp';
    if (gcEl)  gcEl.textContent  = gc + '%';
    
    if (tm !== undefined) {
        const tmElement = document.getElementById(`${primerPrefix}-tm`);
        if (tmElement) tmElement.textContent = tm + '°C';
    }
    
    if (dg5 !== undefined) {
        const dg5Element = document.getElementById(`${primerPrefix}-dg5`);
        if (dg5Element) dg5Element.textContent = dg5 + ' kcal/mol';
    }
    
    if (dg3 !== undefined) {
        const dg3Element = document.getElementById(`${primerPrefix}-dg3`);
        if (dg3Element) dg3Element.textContent = dg3 + ' kcal/mol';
    }
    
    // Run hairpin detection automatically
    const hairpinElement = document.getElementById(`${primerPrefix}-hairpin`);
    if (hairpinElement && sequence) {
        const hairpinResult = detectHairpin(sequence);
        hairpinElement.textContent = hairpinResult.summary;

        if (hairpinResult.summary === "None") {
            hairpinElement.className = 'stat-value good';
        } else if (hairpinResult.strong.length > 0) {
            hairpinElement.className = 'stat-value bad';
        } else {
            hairpinElement.className = 'stat-value warning';
        }
    }
}

// Copy sequence to clipboard
function copySequence(seqId) {
    const seqElement = document.getElementById(seqId);
    const sequence = seqElement.tagName === 'INPUT'
        ? seqElement.value
        : seqElement.getAttribute('data-sequence');
    
    if (!sequence || sequence === '') {
        alert('No sequence to copy');
        return;
    }
    
    navigator.clipboard.writeText(sequence).then(() => {
        const btn = seqElement.parentElement.querySelector('.btn-copy');
        const originalText = btn.textContent;
        btn.textContent = 'Copied!';
        setTimeout(() => {
            btn.textContent = originalText;
        }, 1500);
    });
}

function detectHairpin(sequence) {
    if (!sequence) return { weak: [], mod: [], strong: [], summary: 'None' };  // ← ADD THIS
    sequence = sequence.toUpperCase();

    const results = { weak: [], mod: [], strong: [] };

    // A real hairpin: stem + loop (3–8 nt) + reverse complement of stem
    const MIN_LOOP = 3;
    const MAX_LOOP = 8;

    let i = 0;
    while (i < sequence.length - 2) {
        let found = false;

        // Try stem lengths: strong=5, mod=4, weak=3
        const stemLengths = [
            { len: 5, bucket: 'strong' },
            { len: 4, bucket: 'mod' },
            { len: 3, bucket: 'weak' }
        ];

        for (const { len, bucket } of stemLengths) {
            if (found) break;
            if (i + len + MIN_LOOP + len > sequence.length) continue;

            const stem = sequence.substring(i, i + len);
            const rc   = reverseComplement(stem);
            if (!rc) continue;

            // Check every possible loop size from MIN to MAX
            for (let loop = MIN_LOOP; loop <= MAX_LOOP; loop++) {
                const rcStart = i + len + loop;
                if (rcStart + len > sequence.length) break;

                // RC arm must start exactly at rcStart
                if (sequence.substring(rcStart, rcStart + len) === rc) {
                    results[bucket].push({ seq: stem, pos: i, loopSize: loop });
                    i += len; // skip past stem, no overlaps
                    found = true;
                    break;
                }
            }
        }

        if (!found) i++;
    }

    const total = results.weak.length + results.mod.length + results.strong.length;

    const summary = total === 0
        ? "None"
        : [
            results.strong.length ? `Strong: ${results.strong.length}` : null,
            results.mod.length    ? `Mod: ${results.mod.length}`       : null,
            results.weak.length   ? `Weak: ${results.weak.length}`     : null,
        ].filter(Boolean).join('\n');
    
    return { weak: results.weak, mod: results.mod, strong: results.strong, summary };
}

function detectDimer(sequence1, sequence2) {
    const results = [];
    let i = 0;
 
    while (i < sequence1.length - 2) {
        let found = false;
 
        if (i + 5 <= sequence1.length) {
            const sub5 = sequence1.substring(i, i + 5);
            const rc5  = reverseComplement(sub5);
            if (rc5 && sequence2.includes(rc5)) {
                results.push({ label: `Strong at ${i+1}`, pos: i, len: 5 });
                i += 5;
                found = true;
            }
        }
        if (!found && i + 4 <= sequence1.length) {
            const sub4 = sequence1.substring(i, i + 4);
            const rc4  = reverseComplement(sub4);
            if (rc4 && sequence2.includes(rc4)) {
                results.push({ label: `Mod. at ${i+1}`, pos: i, len: 4 });
                i += 4;
                found = true;
            }
        }
 
        if (!found) i++;
    }
 
    if (results.length === 0) return [{ label: "None", pos: -1, len: 0 }];
    return results;
}

function countAllDimers() {
    const primerNames = Object.keys(designState.outputs).filter(p => p !== 'template');

    // Keys are always nameA:nameB where nameA comes first in primerNames order
    // so we only need one direction per pair here — no need for reverse keys
    const skipPairs = new Set([
        'fip:f2', 'fip:f1c',
        'bip:b2', 'bip:b1c',
        'lf:f2', 'lb:f2', 'f1c:f2', 'b1c:f2', 'f2:b2',
        'lf:b2', 'lb:b2', 'f1c:b2', 'b1c:b2',
        'fip:b2', 'bip:f2'
    ]);

    let total = 0;
    // let hasStrong = false;
    // const pairs = [];

    for (let i = 0; i < primerNames.length; i++) {
        for (let j = i + 1; j < primerNames.length; j++) {
            const nameA = primerNames[i];
            const nameB = primerNames[j];

            const key1 = `${nameA}:${nameB}`;
            const key2 = `${nameB}:${nameA}`;
            if (skipPairs.has(key1) || skipPairs.has(key2)) continue;

            const seqA = designState.outputs[nameA]?.seq;
            const seqB = designState.outputs[nameB]?.seq;
            if (!seqA || !seqB) continue;

            const results = detectDimer(seqA, seqB);
            const realDimers = results.filter(r => r.pos >= 0);

            if (realDimers.length > 0) {
                total += 1;
                // if (realDimers.some(r => r.label.includes("Strong"))) hasStrong = true;
                // pairs.push({ nameA, nameB, dimers: realDimers });
            }
        }
    }

    return { total };
    // return { total, hasStrong, pairs };
}

function updateDimerBox() {
    const el = document.getElementById('dimer-summary');
    if (!el) return;

    const { total } = countAllDimers();

    if (total === 0) {
        el.textContent = 'No primer dimers';
        el.className   = 'stat-value good';
        return;
    }

    el.innerHTML = `<u>${total} dimer pair${total !== 1 ? 's' : ''}</u>`;
}

// Clear form
function clearForm() {
    document.getElementById('mirna1-name').value = '';
    document.getElementById('mirna1-sequence').value = '';
    document.getElementById('mirna2-name').value = '';
    document.getElementById('mirna2-sequence').value = '';
    
    document.getElementById('mirna1-info').classList.remove('visible');
    document.getElementById('mirna2-info').classList.remove('visible');
    document.getElementById('mirna1-error').classList.remove('visible');
    document.getElementById('mirna2-error').classList.remove('visible');
    
    // Reset to default architecture
    selectArchitecture('f2-only');
    document.querySelector('input[value="f2-only"]').checked = true;
}

//Get the reverse complement
function reverseComplement(sequence){
    //Dictionary for bases and their comeplements
    const complement = {
        A: 'T',
        T: 'A',
        C: 'G',
        G: 'C' 
    };
    
     // Normalize input
    sequence = sequence.toUpperCase();

    const bases = sequence.split('');

    //Check if correct
    for (let i = 0; i < bases.length; i++) {
        if (
            bases[i] !== 'A' &&
            bases[i] !== 'T' &&
            bases[i] !== 'C' &&
            bases[i] !== 'G'
        ) {
            console.error("Incorrect sequence for reverse complement");
            return null;
        }
    }

    //.split('') turns into array
    //.reverse reverses the sequence
    //.map 
    //Return to a string
    return bases
        .reverse()
        .map(base => complement[base])
        .join('');
}

function highlightTemplate(rawSequence) {

    // Step 1 — apply primer color spans first, on the clean raw sequence
    const primers = {
        fip: designState.outputs.fip.seq,
        bip: designState.outputs.bip.seq ? reverseComplement(designState.outputs.bip.seq) : '',
        lf:  designState.outputs.lf.seq  ? reverseComplement(designState.outputs.lf.seq)  : '',
        lb:  designState.outputs.lb.seq,
    };

    let annotated = rawSequence;

    const primerList = Object.entries(primers)
        .filter(([_, seq]) => seq)
        .sort((a, b) => b[1].length - a[1].length);

    primerList.forEach(([primer, seq]) => {
        const regex = new RegExp(seq, 'g');
        annotated = annotated.replace(
            regex,
            `<span class="highlight-${primer}" title="${primer.toUpperCase()} primer">${seq}</span>`
        );
    });

    // Step 2 — bold hairpins by replacing raw subseq in the already-spanned string
    const hairpinResult = detectHairpin(rawSequence);
    const allHairpins = [
        ...hairpinResult.strong.map(h => ({ pos: h.pos, len: 5 })),
        ...hairpinResult.mod.map(h =>    ({ pos: h.pos, len: 4 })),
        ...hairpinResult.weak.map(h =>   ({ pos: h.pos, len: 3 })),
    ].sort((a, b) => a.pos - b.pos);

    for (const hp of allHairpins) {
        const subseq = rawSequence.slice(hp.pos, hp.pos + hp.len);
        // Replace only first untagged occurrence to stay accurate
        annotated = annotated.replace(subseq, `<b>${subseq}</b>`);
    }

    // Step 3 — underline dimer regions by finding RC in rawSequence for position,
    // then replacing the raw subseq in the annotated string
    function underlineDimersInTemplate(rawSequence, annotated) {
        // Collect every dimer region (in raw-sequence coordinates) as {start, end}.
        // Use a Set to deduplicate identical coordinate ranges.
        const dimerRanges = new Set();
    
        const allPrimerSeqs = Object.entries(designState.outputs)
            .filter(([name]) => name !== 'template')
            .map(([name, obj]) => ({ name, seq: obj.seq }))
            .filter(p => p.seq);
    
        allPrimerSeqs.forEach(({ name, seq }) => {
            const results = detectDimer(seq, rawSequence);
            results.filter(r => r.pos >= 0).forEach(result => {
                const primerSub = seq.slice(result.pos, result.pos + result.len);
                const rc        = reverseComplement(primerSub);
                if (!rc) return;
    
                // Find ALL occurrences in the template, not just the first.
                let searchFrom = 0;
                while (searchFrom <= rawSequence.length - rc.length) {
                    const hit = rawSequence.indexOf(rc, searchFrom);
                    if (hit < 0) break;
                    dimerRanges.add(`${hit}:${hit + rc.length}:${name}`);
                    searchFrom = hit + 1;   // allow overlapping hits
                }
            });
        });
    
        // Parse the collected ranges into objects and sort by start position.
        const ranges = [...dimerRanges].map(key => {
            const [start, end, name] = key.split(':');
            return { start: +start, end: +end, name };
        }).sort((a, b) => a.start - b.start);
    
        // Inject <u> tags by replacing raw substrings in `annotated`.
        // We replace from right-to-left so earlier string positions stay valid.
        // But since annotated already has <span> tags, we can't use character
        // positions directly — instead we replace each literal raw subseq text
        // carefully, wrapped with a dimer title hint.
        // To avoid double-underlining the same span, track replaced subseqs.
        const replaced = new Set();
        for (const { start, end, name } of ranges) {
            const subseq = rawSequence.slice(start, end);
            const key    = `${start}:${end}`;
            if (replaced.has(key)) continue;
            replaced.add(key);
            // Replace only the FIRST untagged plain occurrence of subseq
            // (i.e., not one already inside a <u> or <b> tag).
            const safeSubseq = subseq.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
            const notInsideTag = new RegExp(`(?<!<[^>]*)${safeSubseq}`, 'g');
            annotated = annotated.replace(notInsideTag, (match) => {
                return `<u title="Dimer with ${name.toUpperCase()}">${match}</u>`;
            });
        }
    
        return annotated;
    }
    annotated = underlineDimersInTemplate(rawSequence, annotated);
    return annotated;
}

//Create template ultramer, disclusing F1, B1, LF, and BF, and spacers
function generateTemplateUltramer(fip, lf, f1c, b1c, lb, bip) {
    const lf_rc = reverseComplement(lf);
    const bip_rc = reverseComplement(bip);
    const f1 = reverseComplement(f1c);

    const rawSequence = fip + lf_rc + 'C' + f1 + 'GT' + b1c + 'G' + lb + 'T' + bip_rc;

    return { rawSequence };
}