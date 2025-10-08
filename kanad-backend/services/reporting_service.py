"""
Reporting Service - LLM-powered report generation using Claude.

This service integrates with Anthropic's Claude API to generate comprehensive,
scientifically accurate reports from computation results.
"""

import anthropic
from typing import Dict, Any, Optional
import logging
import json
import os
from datetime import datetime

logger = logging.getLogger(__name__)


class ReportingService:
    """
    AI-powered report generation service.

    Uses Anthropic's Claude to generate comprehensive scientific reports
    from quantum chemistry computation results.
    """

    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize reporting service.

        Args:
            api_key: Anthropic API key (defaults to env var ANTHROPIC_API_KEY)
        """
        self.api_key = api_key or os.getenv('ANTHROPIC_API_KEY')
        if not self.api_key:
            logger.warning("ANTHROPIC_API_KEY not set. LLM reports will be disabled.")
            self.client = None
        else:
            self.client = anthropic.Anthropic(api_key=self.api_key)
            logger.info("ReportingService initialized with Claude API")

    async def generate_report(
        self,
        molecule: Dict[str, Any],
        results: Dict[str, Any],
        analysis: Dict[str, Any],
        format: str = 'json'
    ) -> Dict[str, Any]:
        """
        Generate comprehensive scientific report.

        Args:
            molecule: Molecule data (formula, name, structure)
            results: Computation results (energy, convergence, etc.)
            analysis: Analysis results (bonds, properties, thermochemistry)
            format: Output format ('json', 'markdown', 'html')

        Returns:
            Generated report dictionary
        """
        if not self.client:
            logger.warning("LLM client not initialized. Returning empty report.")
            return self._generate_fallback_report(molecule, results, analysis)

        logger.info(f"Generating LLM report for {molecule.get('formula', 'Unknown')}")

        try:
            # Construct prompt for Claude
            prompt = self._build_report_prompt(molecule, results, analysis)

            # Call Claude API
            response = self.client.messages.create(
                model="claude-sonnet-4-20250514",  # Latest Sonnet model
                max_tokens=4000,
                temperature=0.3,  # Low temperature for factual, consistent outputs
                system=self._get_system_prompt(),
                messages=[
                    {"role": "user", "content": prompt}
                ]
            )

            # Parse response
            report_text = response.content[0].text

            # Try to parse as JSON
            try:
                report = json.loads(report_text)
            except json.JSONDecodeError:
                # If not JSON, structure manually
                report = {
                    'summary': report_text[:500],
                    'full_text': report_text,
                    'key_findings': self._extract_key_findings(report_text),
                    'interpretation': report_text,
                    'recommendations': []
                }

            # Add metadata
            report['generated_at'] = datetime.utcnow().isoformat()
            report['model'] = 'claude-sonnet-4'
            report['molecule'] = molecule.get('formula')

            logger.info("LLM report generated successfully")
            return report

        except Exception as e:
            logger.error(f"LLM report generation failed: {e}")
            return self._generate_fallback_report(molecule, results, analysis)

    def _get_system_prompt(self) -> str:
        """Get system prompt for Claude."""
        return """You are an expert quantum chemist with deep knowledge of computational chemistry,
quantum mechanics, and molecular properties. Your role is to analyze quantum chemistry computation
results and provide clear, accurate, scientifically rigorous interpretations.

When analyzing results:
1. Be precise with numerical values and units
2. Compare computed values with experimental data when relevant
3. Explain chemical significance of results
4. Identify potential sources of error
5. Suggest follow-up calculations or improvements
6. Use appropriate chemical terminology
7. Be concise but comprehensive

Always structure your response as JSON with these keys:
- summary: Brief 2-3 sentence overview
- key_findings: Array of 3-5 most important results
- interpretation: Detailed scientific interpretation (2-3 paragraphs)
- recommendations: Array of suggested next steps
- comparison_with_experiment: Comparison with known experimental values (if applicable)
- reliability: Assessment of result reliability (high/medium/low) with justification"""

    def _build_report_prompt(
        self,
        molecule: Dict[str, Any],
        results: Dict[str, Any],
        analysis: Dict[str, Any]
    ) -> str:
        """Build prompt for Claude."""
        # Extract key data
        formula = molecule.get('formula', 'Unknown')
        name = molecule.get('name', formula)
        method = results.get('method', 'Unknown')
        energy = results.get('energy', 0.0)
        hf_energy = results.get('hf_energy', 0.0)
        correlation_energy = results.get('correlation_energy', 0.0)
        converged = results.get('converged', False)
        n_iterations = results.get('n_iterations', 0)

        # Build comprehensive prompt
        prompt = f"""Analyze these quantum chemistry computation results:

**Molecule**: {name} ({formula})
**Method**: {method}
**Basis Set**: {molecule.get('basis', 'Unknown')}
**Charge**: {molecule.get('charge', 0)}
**Multiplicity**: {molecule.get('multiplicity', 1)}

**Computational Results**:
- Total Energy: {energy:.6f} Hartree ({energy * 27.2114:.2f} eV)
- HF Energy: {hf_energy:.6f} Hartree
- Correlation Energy: {correlation_energy:.6f} Hartree
- Converged: {converged}
- Iterations: {n_iterations}
"""

        # Add convergence history if available
        if 'convergence_history' in results and results['convergence_history']:
            history = results['convergence_history']
            if len(history) > 0:
                prompt += f"- Initial Energy: {history[0].get('energy', 0):.6f} Ha\n"
                prompt += f"- Final Energy: {history[-1].get('energy', 0):.6f} Ha\n"
                prompt += f"- Energy Change: {abs(history[-1].get('energy', 0) - history[0].get('energy', 0)):.6f} Ha\n"

        # Add analysis results
        if analysis:
            prompt += "\n**Analysis Results**:\n"

            # Energy decomposition
            if 'energy_decomposition' in analysis:
                decomp = analysis['energy_decomposition']
                prompt += f"\nEnergy Decomposition:\n"
                prompt += f"- Kinetic Energy: {decomp.get('kinetic', 0):.6f} Ha\n"
                prompt += f"- Nuclear Attraction: {decomp.get('nuclear_attraction', 0):.6f} Ha\n"
                prompt += f"- Electron Repulsion: {decomp.get('electron_repulsion', 0):.6f} Ha\n"

            # Bond analysis
            if 'bond_analysis' in analysis:
                bonds = analysis['bond_analysis']
                prompt += f"\nBond Analysis:\n"
                prompt += f"- HOMO-LUMO Gap: {bonds.get('homo_lumo_gap', 0):.4f} eV\n"
                if 'bonds' in bonds:
                    prompt += f"- Number of bonds: {len(bonds['bonds'])}\n"

            # Dipole moment
            if 'dipole_moment' in analysis:
                dipole = analysis['dipole_moment']
                prompt += f"\nDipole Moment: {dipole.get('magnitude', 0):.3f} Debye\n"

            # Thermochemistry
            if 'thermochemistry' in analysis:
                thermo = analysis['thermochemistry']
                prompt += f"\nThermochemistry (298K):\n"
                prompt += f"- Enthalpy: {thermo.get('enthalpy', 0):.6f} Ha\n"
                prompt += f"- Entropy: {thermo.get('entropy', 0):.2f} cal/mol·K\n"
                prompt += f"- Gibbs Free Energy: {thermo.get('gibbs_free_energy', 0):.6f} Ha\n"

            # Spectroscopy
            if 'spectroscopy' in analysis:
                spectro = analysis['spectroscopy']
                if 'uv_vis_spectrum' in spectro and spectro['uv_vis_spectrum']:
                    prompt += f"\nUV-Vis Spectrum: {len(spectro['uv_vis_spectrum'])} transitions computed\n"

            # Vibrational analysis
            if 'vibrational' in analysis:
                vib = analysis['vibrational']
                if 'frequencies' in vib and vib['frequencies']:
                    prompt += f"\nVibrational Frequencies: {len(vib['frequencies'])} normal modes\n"
                    prompt += f"- Zero-Point Energy: {vib.get('zero_point_energy', 0):.6f} Ha\n"

            # Excited states
            if 'excited_states' in results:
                states = results['excited_states']
                prompt += f"\nExcited States: {len(states)} states computed\n"
                if states:
                    prompt += f"- First excitation: {states[0].get('excitation_energy', 0):.4f} eV\n"

        prompt += "\n\nProvide a comprehensive scientific analysis in JSON format as specified."

        return prompt

    def _extract_key_findings(self, text: str) -> list:
        """Extract key findings from unstructured text."""
        # Simple extraction - look for bullet points or numbered lists
        findings = []
        for line in text.split('\n'):
            line = line.strip()
            if line.startswith(('-', '*', '•')) or (len(line) > 0 and line[0].isdigit() and '.' in line[:3]):
                findings.append(line.lstrip('-*•0123456789. '))
            if len(findings) >= 5:
                break
        return findings if findings else ["Analysis completed successfully"]

    def _generate_fallback_report(
        self,
        molecule: Dict[str, Any],
        results: Dict[str, Any],
        analysis: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Generate basic report without LLM.

        Used when Claude API is unavailable.
        """
        logger.info("Generating fallback report (no LLM)")

        method = results.get('method', 'Unknown')
        energy = results.get('energy', 0.0)
        formula = molecule.get('formula', 'Unknown')

        summary = f"Computation completed using {method} method on {formula}. "
        summary += f"Final energy: {energy:.6f} Hartree. "
        if results.get('converged'):
            summary += "Calculation converged successfully."
        else:
            summary += "Calculation did not fully converge."

        key_findings = [
            f"Total energy: {energy:.6f} Ha ({energy * 27.2114:.2f} eV)",
            f"Method: {method}",
            f"Convergence: {'Yes' if results.get('converged') else 'No'}"
        ]

        if analysis.get('dipole_moment'):
            dipole = analysis['dipole_moment']['magnitude']
            key_findings.append(f"Dipole moment: {dipole:.3f} D")

        if analysis.get('bond_analysis', {}).get('homo_lumo_gap'):
            gap = analysis['bond_analysis']['homo_lumo_gap']
            key_findings.append(f"HOMO-LUMO gap: {gap:.4f} eV")

        interpretation = f"The {method} calculation on {formula} has been completed. "
        interpretation += f"The computed total energy is {energy:.6f} Hartree. "

        if analysis.get('energy_decomposition'):
            decomp = analysis['energy_decomposition']
            interpretation += f"Energy decomposition shows kinetic energy of {decomp.get('kinetic', 0):.3f} Ha, "
            interpretation += f"nuclear attraction of {decomp.get('nuclear_attraction', 0):.3f} Ha, "
            interpretation += f"and electron repulsion of {decomp.get('electron_repulsion', 0):.3f} Ha. "

        recommendations = [
            "Consider running higher-level method (e.g., CCSD(T)) for more accurate results",
            "Try larger basis set for better energy convergence",
            "Perform frequency analysis to confirm structure is a minimum"
        ]

        return {
            'summary': summary,
            'key_findings': key_findings,
            'interpretation': interpretation,
            'recommendations': recommendations,
            'generated_at': datetime.utcnow().isoformat(),
            'model': 'fallback',
            'molecule': formula
        }

    async def generate_html_report(
        self,
        molecule: Dict[str, Any],
        results: Dict[str, Any],
        analysis: Dict[str, Any]
    ) -> str:
        """
        Generate HTML report.

        Returns:
            HTML string for downloading/viewing
        """
        report = await self.generate_report(molecule, results, analysis)

        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Computation Report - {molecule.get('formula', 'Unknown')}</title>
    <style>
        body {{ font-family: Arial, sans-serif; max-width: 800px; margin: 0 auto; padding: 20px; }}
        h1 {{ color: #2c3e50; }}
        h2 {{ color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 5px; }}
        .summary {{ background: #ecf0f1; padding: 15px; border-left: 4px solid #3498db; margin: 20px 0; }}
        .finding {{ background: #ffffff; padding: 10px; margin: 10px 0; border-left: 3px solid #2ecc71; }}
        .metadata {{ color: #7f8c8d; font-size: 0.9em; }}
        ul {{ line-height: 1.8; }}
    </style>
</head>
<body>
    <h1>Quantum Chemistry Computation Report</h1>
    <p class="metadata">Generated: {report.get('generated_at', 'Unknown')}</p>
    <p class="metadata">Molecule: {molecule.get('name', molecule.get('formula', 'Unknown'))}</p>
    <p class="metadata">Method: {results.get('method', 'Unknown')}</p>

    <h2>Summary</h2>
    <div class="summary">
        <p>{report.get('summary', 'No summary available')}</p>
    </div>

    <h2>Key Findings</h2>
    <ul>
        {''.join([f'<li>{finding}</li>' for finding in report.get('key_findings', [])])}
    </ul>

    <h2>Scientific Interpretation</h2>
    <p>{report.get('interpretation', 'No interpretation available')}</p>

    <h2>Recommendations</h2>
    <ul>
        {''.join([f'<li>{rec}</li>' for rec in report.get('recommendations', [])])}
    </ul>

    <h2>Computational Details</h2>
    <ul>
        <li>Total Energy: {results.get('energy', 0):.6f} Hartree</li>
        <li>HF Energy: {results.get('hf_energy', 0):.6f} Hartree</li>
        <li>Correlation Energy: {results.get('correlation_energy', 0):.6f} Hartree</li>
        <li>Converged: {'Yes' if results.get('converged') else 'No'}</li>
        <li>Iterations: {results.get('n_iterations', 0)}</li>
    </ul>

    <hr>
    <p class="metadata">Report generated by Kanad Quantum Chemistry Platform using Claude AI</p>
</body>
</html>"""

        return html

    async def generate_markdown_report(
        self,
        molecule: Dict[str, Any],
        results: Dict[str, Any],
        analysis: Dict[str, Any]
    ) -> str:
        """Generate Markdown report."""
        report = await self.generate_report(molecule, results, analysis)

        md = f"""# Quantum Chemistry Computation Report

**Generated**: {report.get('generated_at', 'Unknown')}
**Molecule**: {molecule.get('name', molecule.get('formula', 'Unknown'))}
**Method**: {results.get('method', 'Unknown')}

## Summary

{report.get('summary', 'No summary available')}

## Key Findings

"""
        for finding in report.get('key_findings', []):
            md += f"- {finding}\n"

        md += f"""
## Scientific Interpretation

{report.get('interpretation', 'No interpretation available')}

## Recommendations

"""
        for rec in report.get('recommendations', []):
            md += f"- {rec}\n"

        md += f"""
## Computational Details

- **Total Energy**: {results.get('energy', 0):.6f} Hartree
- **HF Energy**: {results.get('hf_energy', 0):.6f} Hartree
- **Correlation Energy**: {results.get('correlation_energy', 0):.6f} Hartree
- **Converged**: {'Yes' if results.get('converged') else 'No'}
- **Iterations**: {results.get('n_iterations', 0)}

---

*Report generated by Kanad Quantum Chemistry Platform using Claude AI*
"""

        return md
