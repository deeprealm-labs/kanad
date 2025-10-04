// API Base URL
const API_BASE = '/api';

// DOM Elements
const atomsContainer = document.getElementById('atoms-container');
const addAtomBtn = document.getElementById('add-atom-btn');
const computeBtn = document.getElementById('compute-btn');
const optimizeBtn = document.getElementById('optimize-btn');
const loading = document.getElementById('loading');
const results = document.getElementById('results');
const error = document.getElementById('error');

// Preset Molecules
const PRESETS = {
    'H2': {
        atoms: [
            { symbol: 'H', position: [0.0, 0.0, 0.0] },
            { symbol: 'H', position: [0.74, 0.0, 0.0] }
        ],
        bondType: 'covalent'
    },
    'LiH': {
        atoms: [
            { symbol: 'Li', position: [0.0, 0.0, 0.0] },
            { symbol: 'H', position: [1.59, 0.0, 0.0] }
        ],
        bondType: 'ionic'
    },
    'H2O': {
        atoms: [
            { symbol: 'O', position: [0.0, 0.0, 0.0] },
            { symbol: 'H', position: [0.96, 0.0, 0.0] },
            { symbol: 'H', position: [-0.24, 0.93, 0.0] }
        ],
        bondType: 'covalent'
    },
    'Na2': {
        atoms: [
            { symbol: 'Na', position: [0.0, 0.0, 0.0] },
            { symbol: 'Na', position: [3.66, 0.0, 0.0] }
        ],
        bondType: 'metallic'
    }
};

// Initialize
document.addEventListener('DOMContentLoaded', () => {
    setupPresetButtons();
    setupAtomManagement();
    setupComputeButtons();
});

// Preset Buttons
function setupPresetButtons() {
    const presetButtons = document.querySelectorAll('.preset-btn');
    presetButtons.forEach(btn => {
        btn.addEventListener('click', () => {
            const presetName = btn.dataset.preset;
            loadPreset(presetName);

            // Visual feedback
            presetButtons.forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
        });
    });
}

function loadPreset(name) {
    const preset = PRESETS[name];
    if (!preset) return;

    // Clear existing atoms
    atomsContainer.innerHTML = '';

    // Add preset atoms
    preset.atoms.forEach(atom => {
        addAtomEntry(atom.symbol, atom.position);
    });

    // Set bond type
    document.getElementById('bond-type').value = preset.bondType;
}

// Atom Management
function setupAtomManagement() {
    addAtomBtn.addEventListener('click', () => {
        addAtomEntry('H', [0.0, 0.0, 0.0]);
    });
}

function addAtomEntry(symbol = 'H', position = [0.0, 0.0, 0.0]) {
    const atomEntry = document.createElement('div');
    atomEntry.className = 'atom-entry';
    atomEntry.innerHTML = `
        <input type="text" class="atom-symbol" placeholder="Symbol" value="${symbol}">
        <input type="number" class="atom-x" placeholder="X" value="${position[0]}" step="0.01">
        <input type="number" class="atom-y" placeholder="Y" value="${position[1]}" step="0.01">
        <input type="number" class="atom-z" placeholder="Z" value="${position[2]}" step="0.01">
        <button class="remove-atom-btn">×</button>
    `;

    const removeBtn = atomEntry.querySelector('.remove-atom-btn');
    removeBtn.addEventListener('click', () => {
        if (atomsContainer.children.length > 1) {
            atomEntry.remove();
        } else {
            showError('At least one atom is required');
        }
    });

    atomsContainer.appendChild(atomEntry);
}

// Compute Buttons
function setupComputeButtons() {
    computeBtn.addEventListener('click', handleCompute);
    optimizeBtn.addEventListener('click', handleOptimize);
}

async function handleCompute() {
    const atoms = getAtoms();
    if (!atoms) return;

    const request = {
        atoms: atoms,
        bond_type: document.getElementById('bond-type').value,
        solver: document.getElementById('solver').value,
        optimize: document.getElementById('optimize').checked,
        strategy: document.getElementById('strategy').value
    };

    showLoading();

    try {
        const response = await fetch(`${API_BASE}/compute`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(request)
        });

        if (!response.ok) {
            const errorData = await response.json();
            throw new Error(errorData.detail || 'Computation failed');
        }

        const data = await response.json();
        displayResults(data);
    } catch (err) {
        showError(err.message);
    }
}

async function handleOptimize() {
    const atoms = getAtoms();
    if (!atoms) return;

    const request = {
        atoms: atoms,
        bond_type: document.getElementById('bond-type').value,
        strategy: document.getElementById('strategy').value
    };

    showLoading();

    try {
        const response = await fetch(`${API_BASE}/optimize`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(request)
        });

        if (!response.ok) {
            const errorData = await response.json();
            throw new Error(errorData.detail || 'Optimization failed');
        }

        const data = await response.json();
        displayOptimizationResults(data);
    } catch (err) {
        showError(err.message);
    }
}

// Get atoms from form
function getAtoms() {
    const atomEntries = atomsContainer.querySelectorAll('.atom-entry');
    const atoms = [];

    for (const entry of atomEntries) {
        const symbol = entry.querySelector('.atom-symbol').value.trim();
        const x = parseFloat(entry.querySelector('.atom-x').value);
        const y = parseFloat(entry.querySelector('.atom-y').value);
        const z = parseFloat(entry.querySelector('.atom-z').value);

        if (!symbol) {
            showError('All atoms must have a chemical symbol');
            return null;
        }

        if (isNaN(x) || isNaN(y) || isNaN(z)) {
            showError('All position coordinates must be valid numbers');
            return null;
        }

        // Proper chemical symbol capitalization: First letter uppercase, rest lowercase
        const properSymbol = symbol.charAt(0).toUpperCase() + symbol.slice(1).toLowerCase();

        atoms.push({
            symbol: properSymbol,
            position: [x, y, z]
        });
    }

    return atoms;
}

// Display Results
function displayResults(data) {
    hideLoading();
    results.classList.remove('hidden');

    // Clear any previous warnings
    const oldWarning = document.querySelector('.warning-message');
    if (oldWarning) oldWarning.remove();

    // Energy
    const energy = data.result.energy;
    document.getElementById('energy-value').textContent = `${energy.toFixed(4)} eV`;
    document.getElementById('method-badge').textContent = data.result.method;

    // System Info
    document.getElementById('orbitals').textContent = data.system.orbitals;
    document.getElementById('electrons').textContent = data.system.electrons;
    document.getElementById('qubits').textContent = data.system.qubits;
    document.getElementById('bond-type-display').textContent = data.system.bond_type;

    // Optimization
    const optimizationCard = document.getElementById('optimization-card');
    if (data.optimization) {
        optimizationCard.classList.remove('hidden');

        const opt = data.optimization;
        const qubitsOriginal = opt.original_qubits || data.system.qubits;
        const qubitsReduced = opt.reduced_qubits || data.system.qubits;
        const reduction = `${qubitsOriginal} → ${qubitsReduced}`;

        document.getElementById('qubit-reduction').textContent = reduction;
        document.getElementById('speedup').textContent = `${opt.speedup_estimate?.toFixed(1) || '1.0'}×`;
        document.getElementById('gate-reduction').textContent = `${opt.gate_reduction_factor?.toFixed(1) || '1.0'}×`;
    } else {
        optimizationCard.classList.add('hidden');
    }

    // Computation Time
    document.getElementById('comp-time').textContent = data.computation_time;

    // Show warning if present (for multi-atom molecules)
    if (data.warning) {
        const warningDiv = document.createElement('div');
        warningDiv.className = 'warning-message';
        warningDiv.textContent = '⚠️ ' + data.warning;
        warningDiv.style.padding = '10px';
        warningDiv.style.background = '#fff3cd';
        warningDiv.style.color = '#856404';
        warningDiv.style.borderRadius = '6px';
        warningDiv.style.marginTop = '10px';
        warningDiv.style.fontSize = '0.9em';
        results.insertBefore(warningDiv, results.firstChild);
    }
}

function displayOptimizationResults(data) {
    hideLoading();
    results.classList.remove('hidden');

    // Show optimization analysis
    const optimizationCard = document.getElementById('optimization-card');
    optimizationCard.classList.remove('hidden');

    const opt = data.optimization;
    const qubitsOriginal = opt.original_qubits || data.system.original_qubits;
    const qubitsReduced = opt.reduced_qubits || data.system.original_qubits;
    const reduction = `${qubitsOriginal} → ${qubitsReduced}`;

    document.getElementById('qubit-reduction').textContent = reduction;
    document.getElementById('speedup').textContent = `${opt.speedup_estimate?.toFixed(1) || '1.0'}×`;
    document.getElementById('gate-reduction').textContent = `${opt.gate_reduction_factor?.toFixed(1) || '1.0'}×`;

    // System info
    document.getElementById('orbitals').textContent = data.system.orbitals;
    document.getElementById('electrons').textContent = data.system.electrons;
    document.getElementById('qubits').textContent = data.system.original_qubits;
    document.getElementById('bond-type-display').textContent = 'Analysis';

    // Hide energy (not computed in optimization-only mode)
    document.getElementById('energy-value').textContent = 'N/A';
    document.getElementById('method-badge').textContent = 'Optimization Analysis';

    document.querySelector('.computation-time').classList.add('hidden');
}

// UI State Management
function showLoading() {
    loading.classList.remove('hidden');
    results.classList.add('hidden');
    error.classList.add('hidden');
}

function hideLoading() {
    loading.classList.add('hidden');
}

function showError(message) {
    hideLoading();
    results.classList.add('hidden');
    error.classList.remove('hidden');
    error.textContent = `Error: ${message}`;
}
