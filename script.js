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
    parameters: {
        targetGC: 50,
        targetTm: 60,
        lengthAdjustment: 0
    }
};

// Normalize sequence parsing - single function used everywhere
function parseSequence(input) {
    // Remove FASTA header(s)
    let cleaned = input.replace(/^>.*$/gm, '');
    
    // Strip whitespace, numbers, and line breaks
    cleaned = cleaned.replace(/[^a-zA-Z]/g, '');
    
    // Convert to uppercase
    cleaned = cleaned.toUpperCase();
    
    // Validate alphabet - only ACGTU allowed
    const validBases = /^[ACGU]+$/;
    if (cleaned && !validBases.test(cleaned)) {
        return { valid: false, sequence: '', error: 'Invalid characters. Only A, C, G, U allowed.' };
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
    
    // Simulate processing (replace with actual Python algorithm call)
    setTimeout(() => {
        // TODO: Call your Python(OR Javascript) algorithm here*************************
        // For now, using mock data*************************
        
        // Update target info
        document.getElementById('result-architecture').textContent = 
            designState.architecture === 'f2-only' ? 'F2 Only' :
            designState.architecture === 'b2-only' ? 'B2 Only' :
            'F2 and B2 (AND-gate)';
        
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
        let fip_sequence = 'CGG​AGA​GGT​CGC​GAT​AGT​CAT' + mirna1Result.sequence;

        updatePrimerOutput('template-seq', 'CGGAGAGGTCGCGATAGTCA...', 150, 55);
        updatePrimerOutput('lf-seq', 'TCACTGATCTGGCCGTAGACCA', 22, 50, 62, -8.5, 'None');
        updatePrimerOutput('lb-seq', 'TGACAGGACATCGGTGACAGT', 21, 52, 61, -7.2, 'None');
        updatePrimerOutput('fip-seq', fip_sequence, fip_sequence.length, calculateGC(fip_sequence), 65, -12.3, '1 weak');
        updatePrimerOutput('bip-seq', 'GATGACAGTGACATCCTGCCTAGGCAGTGTCTTAGCTGGTTGT', 44, 52, 66, -11.8, 'None');
        updatePrimerOutput('f2-seq', mirna1Result.sequence, mirna1Result.sequence.length, calculateGC(mirna1Result.sequence), 58, -6.5);
        updatePrimerOutput('b2-seq', 'TGGCAGTGTCTTAGCTGGTTGT', 22, 50, 59, -7.1);
        updatePrimerOutput('f1c-seq', 'CGGAGAGGTCGCGATAGTCA', 20, 60, 61, -8.2);
        updatePrimerOutput('b1c-seq', 'GATGACAGTGACATCCTGCCT', 21, 52, 60, -7.8);
        
        // Populate customize tab
        document.getElementById('edit-mirna1').value = mirna1Result.sequence;
        if (designState.architecture === 'f2-and-b2') {
            const mirna2Result = parseSequence(document.getElementById('mirna2-sequence').value);
            document.getElementById('edit-mirna2').value = mirna2Result.sequence;
            document.getElementById('edit-mirna2-group').style.display = 'block';
        } else {
            document.getElementById('edit-mirna2-group').style.display = 'none';
        }
        
        // Show results
        setState('results');
    }, 1500);
}

// Update primer output with stable IDs
function updatePrimerOutput(seqId, sequence, length, gc, tm, dg, hairpin) {
    const seqElement = document.getElementById(seqId);
    const primerPrefix = seqId.replace('-seq', '');
    
    // Update sequence
    seqElement.textContent = sequence;
    seqElement.setAttribute('data-sequence', sequence);
    
    // Enable copy button
    const copyBtn = seqElement.parentElement.querySelector('.btn-copy');
    if (copyBtn) copyBtn.disabled = false;
    
    // Update stats
    document.getElementById(`${primerPrefix}-len`).textContent = length + ' bp';
    document.getElementById(`${primerPrefix}-gc`).textContent = gc + '%';
    
    if (tm !== undefined) {
        const tmElement = document.getElementById(`${primerPrefix}-tm`);
        if (tmElement) tmElement.textContent = tm + '°C';
    }
    
    if (dg !== undefined) {
        const dgElement = document.getElementById(`${primerPrefix}-dg`);
        if (dgElement) dgElement.textContent = dg + ' kcal/mol';
    }
    
    if (hairpin !== undefined) {
        const hairpinElement = document.getElementById(`${primerPrefix}-hairpin`);
        if (hairpinElement) {
            hairpinElement.textContent = hairpin;
            hairpinElement.className = 'stat-value ' + 
                (hairpin === 'None' ? 'good' : hairpin.includes('weak') ? 'warning' : 'bad');
        }
    }
}

// Copy sequence to clipboard
function copySequence(seqId) {
    const seqElement = document.getElementById(seqId);
    const sequence = seqElement.getAttribute('data-sequence');
    
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

// Slider updates
function updateGCValue(value) {
    document.getElementById('gc-value').textContent = value + '%';
    designState.parameters.targetGC = parseInt(value);
}

function updateTmValue(value) {
    document.getElementById('tm-value').textContent = value + '°C';
    designState.parameters.targetTm = parseInt(value);
}

function updateLengthValue(value) {
    const sign = value > 0 ? '+' : '';
    document.getElementById('length-value').textContent = sign + value + ' bp';
    designState.parameters.lengthAdjustment = parseInt(value);
}

// Update primers when microRNA changes
function updatePrimers() {
    const newSeq1 = document.getElementById('edit-mirna1').value;
    const result = parseSequence(newSeq1);
    
    if (!result.valid) {
        alert('Invalid sequence: ' + result.error);
        return;
    }
    
    alert('Updating all primers based on new microRNA sequence(s).\n\nThis will recalculate all primer components.');
    // TODO: Call algorithm here to regenerate primers with new sequence****************???
}

// Apply customizations
function applyCustomizations() {
    alert('Applying customizations...\n\nTarget GC: ' + designState.parameters.targetGC + '%\n' +
          'Target Tm: ' + designState.parameters.targetTm + '°C\n' +
          'Length adjustment: ' + designState.parameters.lengthAdjustment + ' bp');
    // TODO: Rerun algorithm with new parameters***********************************************??????
}

// Reset to default
function resetToDefault() {
    document.getElementById('gc-slider').value = 50;
    document.getElementById('tm-slider').value = 60;
    document.getElementById('length-slider').value = 0;
    updateGCValue(50);
    updateTmValue(60);
    updateLengthValue(0);
}