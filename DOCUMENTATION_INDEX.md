# Kanad Documentation Index

This document provides a guide to all documentation and key files in the Kanad project.

## Quick Start Documents

### Start Here (15 minutes)
1. **QUICK_REFERENCE.md** - Project structure, metrics, endpoints, configuration
   - Best for: Quick overview, file locations, API reference
   - Read time: 15 minutes

2. **PROJECT_ARCHITECTURE_OVERVIEW.md** - Comprehensive architecture guide (755 lines)
   - Best for: Understanding all three layers (framework, API, frontend)
   - Read time: 45 minutes

### For Understanding Phase 1 (45 minutes)
3. **PHASE1_ENHANCED_COMPLETE_SUMMARY.md** - Phase 1 achievements, metrics, impact
   - Covers: Hi-VQE, Krylov-SQD, active space, applications, testing
   - Read time: 30 minutes

4. **COMPLETE_HIVQE_STACK_SUMMARY.md** - Hi-VQE implementation details
   - Covers: Architecture, performance metrics, cost analysis
   - Read time: 20 minutes

## Feature Documentation

### Hi-VQE Mode (Measurement Reduction)
- **HIVQE_PROGRESS_SUMMARY.md** - Development progress and achievements
- **HIVQE_QUICK_START_GUIDE.md** - How to use Hi-VQE mode
- **HIVQE_VQE_INTEGRATION_COMPLETE.md** - Integration details
- **HIVQE_ACCURACY_ROOT_CAUSE.md** - Accuracy analysis

### Krylov-SQD Method
- See **COMPLETE_HIVQE_STACK_SUMMARY.md** - Section on Krylov-SQD
- Implemented in: `/home/mk/deeprealm/kanad/kanad/solvers/krylov_sqd_solver.py`

### Active Space Reduction
- **ACTIVE_SPACE_INTEGRATION_COMPLETE.md** - Implementation details
- See configuration options in **QUICK_REFERENCE.md**

### Application Platforms
- **DOMAIN_APPLICATIONS_ARCHITECTURE.md** - Overall application architecture
- **APPLICATIONS_COMPLETE.md** - Completion status and endpoints
- **SESSION_SUMMARY_DOMAIN_APPLICATIONS.md** - Application platform work

## Technical Reference

### Framework Architecture
- **kanad/** - Main framework directory
  - **solvers/** - Quantum solvers (VQE, SQD, Krylov-SQD, excited states)
  - **ansatze/** - Quantum circuits (9 types)
  - **backends/** - Cloud backends (IBM, BlueQubit)
  - **applications/** - Domain platforms (4 types)
  - **analysis/** - Analysis tools (energy, ADME, spectroscopy, etc.)
  - **governance/** - Physics-aware protocols
  - **core/** - Molecules, Hamiltonians, mappers

### API Routes
- **api/main.py** (173 lines) - FastAPI server
- **api/routes/configuration.py** (25 KB) - Framework options
- **api/routes/applications.py** (700 lines) - Domain platforms
- **api/routes/experiments.py** - Experiment management
- **api/services/experiment_service.py** (1355 lines) - Execution engine
- **api/services/application_service.py** - Application service layer

### Frontend Components
- **web/src/components/** (54+ components)
  - **dashboard/** - Dashboard and main UI
  - **simulation/** - Experiment config and results
  - **settings/** - User and backend settings
  - **molecule/** - Molecule creation and viewing
  - **ui/** - Reusable UI elements

## Status & Progress Documents

### Current Phase Status
- **SESSION_STATUS.md** - Latest session status
- **SESSION_STATUS_NOV_4.md** - November 4 session summary
- **REMAINING_WORK_PLAN.md** - Work plan and priorities

### Testing & Validation
- **VALIDATION_RESULTS.md** - Test results and validation
- Tests located in: `/home/mk/deeprealm/kanad/tests/`

## Deployment & Configuration

### Deployment Guides
- **DEPLOYMENT.md** - Deployment instructions
- **IBM_HARDWARE_DEPLOYMENT_GUIDE.md** - IBM Quantum deployment
- **SUCCESSFUL_IBM_DEPLOYMENT.md** - Deployment success report
- **REAL_QISKIT_CIRCUIT_IMPLEMENTATION.md** - Qiskit integration

### Configuration
- **.env** - Environment variables (API keys, database, etc.)
- **.env.example** - Example configuration

## Business & Strategy

### Market Analysis
- **KANAD_COMPETITIVE_STRATEGY.md** - Competitive positioning
- See cost analysis in **QUICK_REFERENCE.md**

### Error Handling
- **ERROR_MITIGATION_STRATEGY.md** - Error mitigation approaches
- **SQD_NOISE_MITIGATION_STRATEGY.md** - SQD-specific noise mitigation
- **HIVQE_ACCURACY_ROOT_CAUSE.md** - Hi-VQE accuracy analysis

## Other Documentation

### Architecture Diagrams
- **TWO_LAYER_ARCHITECTURE.md** - Two-layer system design
- **DOMAIN_APPLICATIONS_ARCHITECTURE.md** - Application architecture

### Backend Integration
- **IBM_BACKEND_INTEGRATION_COMPLETE.md** - IBM backend details
- **IBM_DEPLOYMENT_STATUS.md** - IBM deployment status

### Comparative Analysis
- **SAMPLER_VS_ESTIMATOR_ANALYSIS.md** - Sampler vs Estimator analysis
- **GRADIENT_SELECTION_SUCCESS.md** - Gradient selection analysis

## File Organization

```
/home/mk/deeprealm/kanad/
├── Documentation Files (35+ .md files)
│   ├── QUICK_REFERENCE.md (START HERE)
│   ├── PROJECT_ARCHITECTURE_OVERVIEW.md
│   ├── PHASE1_ENHANCED_COMPLETE_SUMMARY.md
│   ├── COMPLETE_HIVQE_STACK_SUMMARY.md
│   └── ... (30+ more documentation files)
│
├── kanad/                    # Framework (20 modules)
│   ├── solvers/              # 5 solver types
│   ├── ansatze/              # 9 ansatz types
│   ├── backends/             # IBM, BlueQubit
│   ├── applications/         # 4 domains
│   ├── analysis/             # 10+ analysis tools
│   ├── governance/           # 3 protocols
│   ├── core/                 # Molecules, Hamiltonians
│   └── environment/          # Environmental effects
│
├── api/                      # FastAPI backend
│   ├── routes/               # 10+ endpoint modules
│   ├── services/             # Business logic
│   ├── core/                 # Database, config
│   ├── auth/                 # Authentication
│   └── middleware/           # Security, rate limiting
│
├── web/                      # Next.js frontend
│   └── src/
│       ├── components/       # 54+ React components
│       ├── app/              # Pages
│       ├── store/            # Redux state
│       └── lib/              # API client
│
└── tests/                    # Test suite
```

## How to Navigate This Documentation

### If you want to understand...

**The whole project**:
1. Start with QUICK_REFERENCE.md (15 min)
2. Read PROJECT_ARCHITECTURE_OVERVIEW.md (45 min)
3. Explore kanad/ directory structure

**Hi-VQE mode**:
1. Read QUICK_REFERENCE.md section on Hi-VQE
2. Read HIVQE_QUICK_START_GUIDE.md
3. Read HIVQE_PROGRESS_SUMMARY.md
4. Examine kanad/utils/vqe_solver.py

**Phase 1 achievements**:
1. Read PHASE1_ENHANCED_COMPLETE_SUMMARY.md
2. Read COMPLETE_HIVQE_STACK_SUMMARY.md
3. Review test results in tests/ directory

**API endpoints**:
1. See QUICK_REFERENCE.md API Endpoints section
2. Read api/routes/configuration.py (options endpoint)
3. Read api/routes/applications.py (domain platforms)

**Application platforms**:
1. Read DOMAIN_APPLICATIONS_ARCHITECTURE.md
2. Explore kanad/applications/ directory
3. See api/services/application_service.py

**Frontend architecture**:
1. See PROJECT_ARCHITECTURE_OVERVIEW.md Part 3
2. Explore web/src/components/ directory
3. Check web/src/lib/api.ts for API integration

**Cloud deployment**:
1. Read DEPLOYMENT.md
2. Read IBM_HARDWARE_DEPLOYMENT_GUIDE.md
3. See api/.env.example for configuration

## Key Metrics at a Glance

| Metric | Value |
|--------|-------|
| Framework Utilization | 85% (Phase 1 complete) |
| Test Pass Rate | 96.5% (441/457) |
| Hi-VQE Cost Savings | 99.98% ($15,000 → $3) |
| Hi-VQE Measurement Reduction | 1000x |
| Krylov-SQD Efficiency Gain | 10-20x |
| Active Space Qubit Reduction | 17% |
| API Endpoints | 40+ |
| UI Components | 54+ |
| Code Base | ~50,000 lines |

## Recent Changes Summary

- **Phase 1 Complete**: Hi-VQE, Krylov-SQD, active space, applications
- **Modified Files**: 25+ (framework, API, frontend)
- **Deleted Files**: 6 (Azure configs)
- **Untracked Files**: 40+ (tests, results)
- **Status**: Production-ready, awaiting Phase 2

## Next Steps

1. Review uncommitted Phase 1 changes
2. Commit Phase 1 work
3. Begin Phase 2: Environmental effects integration
4. Add multi-atom Krylov-SQD support
5. Enhance ADAPT-VQE integration

## Contact & Support

For questions about specific components:
- Framework: See kanad/ module docstrings
- API: See api/routes/ and api/services/ docstrings
- Frontend: See web/src/components/ JSDoc comments
- Tests: See tests/ for examples

---

**Last Updated**: November 6, 2025
**Status**: Phase 1 Complete, Production-Ready
**Next Phase**: Environmental Effects Integration
