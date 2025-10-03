# Kanad Test Suite

Comprehensive testing for the Kanad quantum chemistry framework.

## Directory Structure

```
tests/
├── unit/              # Unit tests for individual components
├── integration/       # Integration tests across modules
└── validation/        # Scientific validation and benchmarks
```

## Running Tests

### All Unit Tests
```bash
pytest tests/unit/ -v
```

### Specific Test Module
```bash
pytest tests/unit/test_bonds.py -v
```

### Integration Tests
```bash
pytest tests/integration/ -v
```

### Validation Suite
```bash
python tests/validation/run_all_validations.py
```

## Test Coverage

- **Unit Tests**: 276 passing tests covering all core modules
- **Integration Tests**: Cross-module functionality
- **Validation**: 24 scientific validation scripts
  - Molecular benchmarks (H₂, NaCl, H₂O)
  - Metallic systems (1D/2D/3D lattices)
  - Thermodynamic properties
  - VQE integration
  - Bond comparisons

## CI/CD

Tests are automatically run on:
- Pull requests
- Main branch commits
- Nightly builds (validation suite)
