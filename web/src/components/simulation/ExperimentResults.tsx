"use client";

import { useState } from "react";
import { LineChart, Line, BarChart, Bar, ScatterChart, Scatter, CartesianGrid, XAxis, YAxis, Tooltip, ResponsiveContainer, Legend, Cell } from "recharts";
import { Activity, Atom, Zap, Layers, TrendingDown, TrendingUp } from "lucide-react";

interface ExperimentResultsProps {
  results: any;
  experimentConfig: any;
}

export default function ExperimentResults({ results, experimentConfig }: ExperimentResultsProps) {
  const [activeTab, setActiveTab] = useState<"energy" | "analysis" | "application">("energy");

  // Determine application domain from config
  const applicationDomain = experimentConfig?.applicationDomain || experimentConfig?.application_domain || null;

  console.log("ExperimentResults - applicationDomain:", applicationDomain);
  console.log("ExperimentResults - results:", results);

  return (
    <div className="h-full flex flex-col bg-background overflow-hidden">
      {/* Tabs Navigation */}
      <div className="flex gap-4 px-6 py-3 border-b border-border flex-shrink-0">
        <button
          onClick={() => setActiveTab("energy")}
          className={`px-4 py-2 text-sm font-quando font-semibold rounded-lg transition ${
            activeTab === "energy"
              ? "bg-brand-orange text-white"
              : "text-muted-foreground hover:text-foreground hover:bg-accent"
          }`}
        >
          <div className="flex items-center gap-2">
            <Zap className="w-4 h-4" />
            Energy & Convergence
          </div>
        </button>
        <button
          onClick={() => setActiveTab("analysis")}
          className={`px-4 py-2 text-sm font-quando font-semibold rounded-lg transition ${
            activeTab === "analysis"
              ? "bg-brand-orange text-white"
              : "text-muted-foreground hover:text-foreground hover:bg-accent"
          }`}
        >
          <div className="flex items-center gap-2">
            <Activity className="w-4 h-4" />
            Molecular Analysis
          </div>
        </button>
        {applicationDomain && (
          <button
            onClick={() => setActiveTab("application")}
            className={`px-4 py-2 text-sm font-quando font-semibold rounded-lg transition ${
              activeTab === "application"
                ? "bg-brand-orange text-white"
                : "text-muted-foreground hover:text-foreground hover:bg-accent"
            }`}
          >
            <div className="flex items-center gap-2">
              <Layers className="w-4 h-4" />
              {applicationDomain === "drug_discovery" && "Drug Discovery"}
              {applicationDomain === "materials_science" && "Materials Science"}
              {applicationDomain === "catalysis" && "Catalysis"}
              {applicationDomain === "energy_storage" && "Energy Storage"}
            </div>
          </button>
        )}
      </div>

      {/* Tab Content */}
      <div className="flex-1 overflow-y-auto p-6 min-h-0">
        {activeTab === "energy" && <EnergyTab results={results} />}
        {activeTab === "analysis" && <AnalysisTab results={results} />}
        {activeTab === "application" && applicationDomain && (
          <ApplicationTab results={results} domain={applicationDomain} />
        )}
      </div>
    </div>
  );
}

// Energy & Convergence Tab
function EnergyTab({ results }: { results: any }) {
  const convergenceData = results.convergence_history ||
    (results.energy_history ? results.energy_history.map((e: number, i: number) => ({ iteration: i + 1, energy: e })) : []);

  return (
    <div className="grid grid-cols-2 gap-4 h-full max-h-[calc(100vh-16rem)]">
      {/* Energy Summary */}
      <div className="bg-card border border-border rounded-lg p-4 overflow-y-auto">
        <h3 className="text-base font-quando font-bold mb-4">Energy Summary</h3>
        <div className="space-y-3">
          <div className="bg-gradient-to-r from-brand-orange/10 to-brand-orange/5 rounded-lg p-4 border border-brand-orange/20">
            <div className="text-sm text-muted-foreground mb-1">Ground State Energy</div>
            <div className="text-3xl font-quando font-bold text-brand-orange">
              {results.energy?.toFixed(6) || "N/A"} <span className="text-lg">Ha</span>
            </div>
          </div>

          {results.hf_energy && (
            <div className="bg-muted rounded-lg p-4">
              <div className="text-sm text-muted-foreground mb-1">Hartree-Fock Energy</div>
              <div className="text-2xl font-quanto font-bold">
                {results.hf_energy.toFixed(6)} <span className="text-base">Ha</span>
              </div>
            </div>
          )}

          {results.correlation_energy !== undefined && (
            <div className="bg-muted rounded-lg p-4">
              <div className="text-sm text-muted-foreground mb-1">Correlation Energy</div>
              <div className="text-2xl font-quanto font-bold text-purple-600 dark:text-purple-400">
                {results.correlation_energy.toFixed(6)} <span className="text-base">Ha</span>
              </div>
            </div>
          )}

          {results.nuclear_repulsion !== undefined && (
            <div className="bg-muted rounded-lg p-4">
              <div className="text-sm text-muted-foreground mb-1">Nuclear Repulsion</div>
              <div className="text-2xl font-quanto font-bold">
                {results.nuclear_repulsion.toFixed(6)} <span className="text-base">Ha</span>
              </div>
            </div>
          )}

          <div className="grid grid-cols-2 gap-3 mt-4">
            <div className="bg-muted rounded-lg p-3">
              <div className="text-xs text-muted-foreground mb-1">Iterations</div>
              <div className="text-xl font-quanto font-bold">{results.iterations || "N/A"}</div>
            </div>
            <div className="bg-muted rounded-lg p-3">
              <div className="text-xs text-muted-foreground mb-1">Converged</div>
              <div className="text-xl font-quanto font-bold">
                {results.converged ? "✓ Yes" : "✗ No"}
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* Convergence Graph */}
      <div className="bg-card border border-border rounded-lg p-4 overflow-hidden flex flex-col">
        <h3 className="text-base font-quando font-bold mb-3">Convergence History</h3>
        <div className="flex-1 min-h-0">
          <ResponsiveContainer width="100%" height="100%">
          <LineChart data={convergenceData}>
            <CartesianGrid strokeDasharray="3 3" stroke="hsl(var(--border))" />
            <XAxis
              dataKey="iteration"
              label={{ value: "Iteration", position: "insideBottom", offset: -5 }}
              tick={{ fontSize: 11 }}
              stroke="hsl(var(--muted-foreground))"
            />
            <YAxis
              label={{ value: "Energy (Ha)", angle: -90, position: "insideLeft" }}
              tick={{ fontSize: 11 }}
              stroke="hsl(var(--muted-foreground))"
              domain={["auto", "auto"]}
            />
            <Tooltip
              contentStyle={{
                backgroundColor: "hsl(var(--card))",
                border: "1px solid hsl(var(--border))",
                borderRadius: "8px",
                fontSize: "12px",
              }}
              formatter={(value: any) => [value.toFixed(8) + " Ha", "Energy"]}
            />
            <Line
              type="monotone"
              dataKey="energy"
              stroke="#ea580c"
              strokeWidth={3}
              dot={{ fill: "#ea580c", r: 3 }}
              activeDot={{ r: 6 }}
            />
          </LineChart>
        </ResponsiveContainer>
        </div>
      </div>
    </div>
  );
}

// Molecular Analysis Tab
function AnalysisTab({ results }: { results: any }) {
  const analysis = results.analysis || {};
  const hasOrbitalEnergies = results.orbital_energies && Array.isArray(results.orbital_energies) && results.orbital_energies.length > 0;
  const hasDipole = results.dipole && Array.isArray(results.dipole);

  return (
    <div className="grid grid-cols-2 gap-4 h-full max-h-[calc(100vh-16rem)]">
      {/* Molecular Properties */}
      <div className="bg-card border border-border rounded-lg p-4 overflow-y-auto">
        <h3 className="text-base font-quando font-bold mb-4 flex items-center gap-2">
          <Atom className="w-4 h-4" />
          Molecular Properties
        </h3>
        <div className="space-y-3">
          {results.n_atoms && (
            <div className="bg-muted rounded-lg p-4">
              <div className="text-sm text-muted-foreground mb-1">Number of Atoms</div>
              <div className="text-2xl font-quanto font-bold">{results.n_atoms}</div>
            </div>
          )}

          {results.n_electrons && (
            <div className="bg-muted rounded-lg p-4">
              <div className="text-sm text-muted-foreground mb-1">Number of Electrons</div>
              <div className="text-2xl font-quanto font-bold">{results.n_electrons}</div>
            </div>
          )}

          {results.charge !== undefined && (
            <div className="bg-muted rounded-lg p-4">
              <div className="text-sm text-muted-foreground mb-1">Charge</div>
              <div className="text-2xl font-quanto font-bold">{results.charge}</div>
            </div>
          )}

          {results.multiplicity && (
            <div className="bg-muted rounded-lg p-4">
              <div className="text-sm text-muted-foreground mb-1">Multiplicity</div>
              <div className="text-2xl font-quanto font-bold">{results.multiplicity}</div>
            </div>
          )}

          {hasDipole && (
            <div className="bg-muted rounded-lg p-4">
              <div className="text-sm text-muted-foreground mb-1">Dipole Moment</div>
              <div className="text-xl font-quanto font-bold">
                {Math.sqrt(results.dipole.reduce((sum: number, d: number) => sum + d*d, 0)).toFixed(4)} D
              </div>
              <div className="text-xs text-muted-foreground mt-2 font-mono">
                ({results.dipole[0].toFixed(4)}, {results.dipole[1].toFixed(4)}, {results.dipole[2].toFixed(4)})
              </div>
            </div>
          )}
        </div>
      </div>

      {/* Orbital Energies */}
      {hasOrbitalEnergies && (
        <div className="bg-card border border-border rounded-lg p-4 overflow-y-auto">
          <h3 className="text-base font-quando font-bold mb-4">Orbital Energies</h3>
          <div className="space-y-3">
            <div className="grid grid-cols-2 gap-4">
              <div className="bg-gradient-to-br from-green-50 to-green-100 dark:from-green-900/30 dark:to-green-900/20 border border-green-200 dark:border-green-800 rounded-lg p-4">
                <div className="text-green-700 dark:text-green-300 text-sm mb-2 font-semibold">HOMO</div>
                <div className="text-2xl font-quanto font-bold text-green-900 dark:text-green-100">
                  {results.orbital_energies[results.orbital_energies.length - 1].toFixed(4)}
                </div>
                <div className="text-xs text-green-600 dark:text-green-400 mt-1">eV</div>
              </div>

              {results.orbital_energies.length > 1 && (
                <div className="bg-gradient-to-br from-red-50 to-red-100 dark:from-red-900/30 dark:to-red-900/20 border border-red-200 dark:border-red-800 rounded-lg p-4">
                  <div className="text-red-700 dark:text-red-300 text-sm mb-2 font-semibold">LUMO</div>
                  <div className="text-2xl font-quanto font-bold text-red-900 dark:text-red-100">
                    {results.orbital_energies[0].toFixed(4)}
                  </div>
                  <div className="text-xs text-red-600 dark:text-red-400 mt-1">eV</div>
                </div>
              )}
            </div>

            {results.orbital_energies.length > 1 && (
              <div className="bg-gradient-to-r from-blue-50 to-purple-50 dark:from-blue-900/20 dark:to-purple-900/20 border border-blue-200 dark:border-blue-800 rounded-lg p-4">
                <div className="text-blue-700 dark:text-blue-300 text-sm mb-2 font-semibold">HOMO-LUMO Gap</div>
                <div className="text-3xl font-quanto font-bold text-blue-900 dark:text-blue-100">
                  {Math.abs(results.orbital_energies[0] - results.orbital_energies[results.orbital_energies.length - 1]).toFixed(4)}
                </div>
                <div className="text-xs text-blue-600 dark:text-blue-400 mt-1">eV</div>
              </div>
            )}

            {/* Orbital Energy Diagram */}
            <div className="bg-muted rounded-lg p-4 mt-6">
              <div className="text-sm text-muted-foreground mb-4">Molecular Orbital Diagram</div>
              <ResponsiveContainer width="100%" height={200}>
                <BarChart
                  data={results.orbital_energies.map((e: number, i: number) => ({
                    index: i + 1,
                    energy: e,
                    type: i >= results.orbital_energies.length - (results.n_electrons / 2) ? "occupied" : "virtual"
                  }))}
                  layout="horizontal"
                >
                  <CartesianGrid strokeDasharray="3 3" stroke="hsl(var(--border))" />
                  <XAxis type="number" dataKey="energy" label={{ value: "Energy (eV)", position: "insideBottom", offset: -5 }} />
                  <YAxis type="category" dataKey="index" label={{ value: "Orbital", angle: -90, position: "insideLeft" }} />
                  <Tooltip />
                  <Bar dataKey="energy">
                    {results.orbital_energies.map((e: number, i: number) => (
                      <Cell key={i} fill={i >= results.orbital_energies.length - (results.n_electrons / 2) ? "#22c55e" : "#ef4444"} />
                    ))}
                  </Bar>
                </BarChart>
              </ResponsiveContainer>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}

// Application-Specific Tab
function ApplicationTab({ results, domain }: { results: any; domain: string }) {
  // Get application-specific analysis from backend
  const applicationData = results.application?.analysis || {};

  // Fallback to advanced_analysis for backward compatibility
  const advancedAnalysis = results.advanced_analysis?.results || {};

  // Merge both sources (application data takes precedence)
  const analysis = { ...advancedAnalysis, ...applicationData };

  console.log("ApplicationTab - domain:", domain);
  console.log("ApplicationTab - applicationData:", applicationData);
  console.log("ApplicationTab - analysis:", analysis);

  if (domain === "drug_discovery" || domain === "drug-discovery") {
    return <DrugDiscoveryResults results={results} analysis={analysis} applicationData={applicationData} />;
  } else if (domain === "materials_science" || domain === "materials-science") {
    return <MaterialsScienceResults results={results} analysis={analysis} applicationData={applicationData} />;
  } else if (domain === "catalysis") {
    return <CatalysisResults results={results} analysis={analysis} applicationData={applicationData} />;
  } else if (domain === "energy_storage" || domain === "energy-storage") {
    return <EnergyStorageResults results={results} analysis={analysis} applicationData={applicationData} />;
  }

  return (
    <div className="flex items-center justify-center h-full">
      <div className="text-center">
        <div className="text-muted-foreground mb-2">Application-specific results not available</div>
        <div className="text-sm text-muted-foreground">Domain: {domain}</div>
      </div>
    </div>
  );
}

// Drug Discovery Results
function DrugDiscoveryResults({ results, analysis, applicationData }: { results: any; analysis: any; applicationData: any }) {
  // Extract ADME properties from backend response
  const adme = applicationData?.adme_properties || analysis?.adme_properties || {};
  const quantum = applicationData?.quantum_properties || {};

  // Calculate HOMO-LUMO gap
  const homoLumoGap = quantum.homo_lumo_gap ||
    (results.orbital_energies && results.orbital_energies.length > 1
      ? Math.abs(results.orbital_energies[0] - results.orbital_energies[results.orbital_energies.length - 1])
      : null);

  // Lipinski checks
  const mw = adme.molecular_weight || 0;
  const logP = adme.logP || 0;
  const hbd = adme.hbd || 0;
  const hba = adme.hba || 0;
  const violations = adme.lipinski_violations || 0;
  const drugScore = adme.druglikeness_score || 0;

  return (
    <div className="grid grid-cols-2 gap-6">
      <div className="bg-card border border-border rounded-lg p-6">
        <h3 className="text-lg font-quando font-bold mb-6">Drug-like Properties</h3>
        <div className="space-y-4">
          <div className="bg-gradient-to-r from-blue-50 to-blue-100 dark:from-blue-900/30 dark:to-blue-900/20 border border-blue-200 dark:border-blue-800 rounded-lg p-4">
            <div className="text-blue-700 dark:text-blue-300 text-sm mb-1">Druglikeness Score</div>
            <div className="text-3xl font-quanto font-bold text-blue-900 dark:text-blue-100">
              {drugScore > 0 ? drugScore.toFixed(2) : "N/A"}
            </div>
            <div className="text-xs text-blue-600 dark:text-blue-400 mt-1">
              {drugScore > 0.8 ? "Excellent" : drugScore > 0.6 ? "Good" : drugScore > 0 ? "Fair" : ""}
            </div>
          </div>

          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">Lipophilicity (LogP)</div>
            <div className="text-2xl font-quanto font-bold">
              {logP > 0 ? logP.toFixed(2) : "N/A"}
            </div>
            <div className="text-xs text-muted-foreground mt-1">
              {logP > 5 ? "⚠️ High" : logP > 0 ? "✓ Optimal" : ""}
            </div>
          </div>

          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">Molecular Weight</div>
            <div className="text-2xl font-quanto font-bold">
              {mw > 0 ? mw.toFixed(2) : "N/A"} {mw > 0 && <span className="text-base">g/mol</span>}
            </div>
            <div className="text-xs text-muted-foreground mt-1">
              {mw > 500 ? "⚠️ High" : mw > 0 ? "✓ Good" : ""}
            </div>
          </div>

          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">HOMO-LUMO Gap</div>
            <div className="text-2xl font-quanto font-bold">
              {homoLumoGap ? homoLumoGap.toFixed(4) : "N/A"} {homoLumoGap && <span className="text-base">eV</span>}
            </div>
            <div className="text-xs text-muted-foreground mt-1">
              {homoLumoGap && (homoLumoGap > 3.0 ? "High stability" : "Moderate reactivity")}
            </div>
          </div>
        </div>
      </div>

      <div className="bg-card border border-border rounded-lg p-6">
        <h3 className="text-lg font-quando font-bold mb-6">Pharmacokinetic Properties</h3>
        <div className="space-y-4">
          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-3">Lipinski's Rule of Five</div>
            <div className="grid grid-cols-2 gap-3">
              <div className="text-xs">
                <div className="text-muted-foreground">MW &lt; 500</div>
                <div className={`font-semibold ${mw <= 500 ? "text-green-600" : "text-red-600"}`}>
                  {mw <= 500 ? "✓ Pass" : "✗ Fail"}
                </div>
              </div>
              <div className="text-xs">
                <div className="text-muted-foreground">LogP &lt; 5</div>
                <div className={`font-semibold ${logP <= 5 ? "text-green-600" : "text-red-600"}`}>
                  {logP <= 5 ? "✓ Pass" : "✗ Fail"}
                </div>
              </div>
              <div className="text-xs">
                <div className="text-muted-foreground">H-bond donors &lt; 5</div>
                <div className={`font-semibold ${hbd <= 5 ? "text-green-600" : "text-red-600"}`}>
                  {hbd <= 5 ? "✓ Pass" : "✗ Fail"}
                </div>
              </div>
              <div className="text-xs">
                <div className="text-muted-foreground">H-bond acceptors &lt; 10</div>
                <div className={`font-semibold ${hba <= 10 ? "text-green-600" : "text-red-600"}`}>
                  {hba <= 10 ? "✓ Pass" : "✗ Fail"}
                </div>
              </div>
            </div>
            <div className="mt-3 text-xs text-center">
              <span className="font-semibold">{violations}</span> violations
            </div>
          </div>

          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">TPSA</div>
            <div className="text-2xl font-quanto font-bold">
              {adme.tpsa ? adme.tpsa.toFixed(2) : "N/A"} {adme.tpsa && <span className="text-base">Ų</span>}
            </div>
            <div className="text-xs text-muted-foreground mt-1">
              {adme.tpsa && (adme.tpsa < 140 ? "Good permeability" : "Limited permeability")}
            </div>
          </div>

          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">Rotatable Bonds</div>
            <div className="text-2xl font-quanto font-bold">
              {adme.n_rotatable_bonds || "N/A"}
            </div>
            <div className="text-xs text-muted-foreground mt-1">
              {adme.n_rotatable_bonds && (adme.n_rotatable_bonds < 10 ? "Good flexibility" : "High flexibility")}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}

// Materials Science Results
function MaterialsScienceResults({ results, analysis, applicationData }: { results: any; analysis: any; applicationData: any }) {
  // Extract materials analysis from backend
  const bandgapData = applicationData?.bandgap || {};
  const opticalData = applicationData?.optical_properties || {};
  const ledData = applicationData?.led_properties || {};

  // Band gap value
  const bandgapValue = bandgapData.value ||
    (results.orbital_energies && results.orbital_energies.length > 1
      ? Math.abs(results.orbital_energies[0] - results.orbital_energies[results.orbital_energies.length - 1])
      : null);

  // Determine semiconductor type
  const semiconductorType = bandgapData.type ||
    (bandgapValue && bandgapValue > 2.0 ? "Semiconductor" : bandgapValue && bandgapValue > 0 ? "Conductor" : "Unknown");

  // LED color
  const ledColor = ledData.emission_color || null;
  const ledWavelength = ledData.emission_wavelength || null;
  const ledRgb = ledData.rgb || null;

  return (
    <div className="grid grid-cols-2 gap-6">
      <div className="bg-card border border-border rounded-lg p-6">
        <h3 className="text-lg font-quando font-bold mb-6">Electronic Properties</h3>
        <div className="space-y-4">
          <div className="bg-gradient-to-r from-purple-50 to-purple-100 dark:from-purple-900/30 dark:to-purple-900/20 border border-purple-200 dark:border-purple-800 rounded-lg p-4">
            <div className="text-purple-700 dark:text-purple-300 text-sm mb-1">Band Gap</div>
            <div className="text-3xl font-quanto font-bold text-purple-900 dark:text-purple-100">
              {bandgapValue ? bandgapValue.toFixed(4) : "N/A"}{" "}
              {bandgapValue && <span className="text-lg">eV</span>}
            </div>
            <div className="text-xs text-purple-600 dark:text-purple-400 mt-1">
              {semiconductorType}
            </div>
          </div>

          {bandgapData.suitable_for_solar !== undefined && (
            <div className="bg-muted rounded-lg p-4">
              <div className="text-sm text-muted-foreground mb-1">Solar Cell Suitability</div>
              <div className={`text-2xl font-quanto font-bold ${bandgapData.suitable_for_solar ? "text-green-600" : "text-yellow-600"}`}>
                {bandgapData.suitable_for_solar ? "✓ Suitable" : "Not Optimal"}
              </div>
              <div className="text-xs text-muted-foreground mt-1">
                Ideal range: 1.0-1.8 eV
              </div>
            </div>
          )}

          {bandgapData.suitable_for_led !== undefined && (
            <div className="bg-muted rounded-lg p-4">
              <div className="text-sm text-muted-foreground mb-1">LED Suitability</div>
              <div className={`text-2xl font-quanto font-bold ${bandgapData.suitable_for_led ? "text-green-600" : "text-yellow-600"}`}>
                {bandgapData.suitable_for_led ? "✓ Suitable" : "Not Optimal"}
              </div>
              <div className="text-xs text-muted-foreground mt-1">
                Ideal range: 1.8-3.5 eV
              </div>
            </div>
          )}
        </div>
      </div>

      <div className="bg-card border border-border rounded-lg p-6">
        <h3 className="text-lg font-quando font-bold mb-6">Optical Properties</h3>
        <div className="space-y-4">
          {opticalData.refractive_index && (
            <div className="bg-muted rounded-lg p-4">
              <div className="text-sm text-muted-foreground mb-1">Refractive Index</div>
              <div className="text-2xl font-quanto font-bold">
                {opticalData.refractive_index.toFixed(4)}
              </div>
            </div>
          )}

          {ledColor && (
            <div className="bg-muted rounded-lg p-4">
              <div className="text-sm text-muted-foreground mb-3">LED Emission Color</div>
              <div className="flex items-center gap-3">
                {ledRgb && (
                  <div
                    className="w-12 h-12 rounded-lg border-2 border-border"
                    style={{ backgroundColor: `rgb(${ledRgb[0]}, ${ledRgb[1]}, ${ledRgb[2]})` }}
                  />
                )}
                <div>
                  <div className="text-2xl font-quanto font-bold capitalize">
                    {ledColor}
                  </div>
                  {ledWavelength && (
                    <div className="text-xs text-muted-foreground">
                      {ledWavelength.toFixed(1)} nm
                    </div>
                  )}
                </div>
              </div>
            </div>
          )}

          {opticalData.absorption_spectrum && (
            <div className="bg-muted rounded-lg p-4">
              <div className="text-sm text-muted-foreground mb-1">Absorption Edge</div>
              <div className="text-2xl font-quanto font-bold">
                {opticalData.absorption_spectrum.absorption_edge?.toFixed(1) || "N/A"} nm
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}

// Catalysis Results
function CatalysisResults({ results, analysis, applicationData }: { results: any; analysis: any; applicationData: any }) {
  return (
    <div className="grid grid-cols-2 gap-6">
      <div className="bg-card border border-border rounded-lg p-6">
        <h3 className="text-lg font-quando font-bold mb-6">Catalytic Activity</h3>
        <div className="space-y-4">
          <div className="bg-gradient-to-r from-orange-50 to-orange-100 dark:from-orange-900/30 dark:to-orange-900/20 border border-orange-200 dark:border-orange-800 rounded-lg p-4">
            <div className="text-orange-700 dark:text-orange-300 text-sm mb-1">Activation Energy</div>
            <div className="text-3xl font-quanto font-bold text-orange-900 dark:text-orange-100">
              {analysis.activation_energy?.toFixed(4) || "N/A"} <span className="text-lg">eV</span>
            </div>
          </div>

          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">Reaction Rate Constant</div>
            <div className="text-2xl font-quanto font-bold">
              {analysis.rate_constant?.toExponential(2) || "N/A"} <span className="text-base">s⁻¹</span>
            </div>
          </div>

          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">Turnover Frequency</div>
            <div className="text-2xl font-quanto font-bold">
              {analysis.tof?.toFixed(2) || "N/A"} <span className="text-base">h⁻¹</span>
            </div>
          </div>
        </div>
      </div>

      <div className="bg-card border border-border rounded-lg p-6">
        <h3 className="text-lg font-quando font-bold mb-6">Surface Properties</h3>
        <div className="space-y-4">
          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">Adsorption Energy</div>
            <div className="text-2xl font-quanto font-bold">
              {analysis.adsorption_energy?.toFixed(4) || "N/A"} <span className="text-base">eV</span>
            </div>
          </div>

          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">Binding Site Type</div>
            <div className="text-2xl font-quanto font-bold text-blue-600 dark:text-blue-400">
              {analysis.binding_site || "Active Site"}
            </div>
          </div>

          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">Selectivity</div>
            <div className="text-2xl font-quanto font-bold text-green-600 dark:text-green-400">
              {analysis.selectivity || "High"}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}

// Energy Storage Results
function EnergyStorageResults({ results, analysis, applicationData }: { results: any; analysis: any; applicationData: any }) {
  return (
    <div className="grid grid-cols-2 gap-6">
      <div className="bg-card border border-border rounded-lg p-6">
        <h3 className="text-lg font-quando font-bold mb-6">Battery Performance</h3>
        <div className="space-y-4">
          <div className="bg-gradient-to-r from-green-50 to-green-100 dark:from-green-900/30 dark:to-green-900/20 border border-green-200 dark:border-green-800 rounded-lg p-4">
            <div className="text-green-700 dark:text-green-300 text-sm mb-1">Voltage</div>
            <div className="text-3xl font-quanto font-bold text-green-900 dark:text-green-100">
              {analysis.voltage?.toFixed(2) || "N/A"} <span className="text-lg">V</span>
            </div>
          </div>

          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">Specific Capacity</div>
            <div className="text-2xl font-quanto font-bold">
              {analysis.specific_capacity?.toFixed(2) || "N/A"} <span className="text-base">mAh/g</span>
            </div>
          </div>

          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">Energy Density</div>
            <div className="text-2xl font-quanto font-bold">
              {analysis.energy_density?.toFixed(2) || "N/A"} <span className="text-base">Wh/kg</span>
            </div>
          </div>
        </div>
      </div>

      <div className="bg-card border border-border rounded-lg p-6">
        <h3 className="text-lg font-quanto font-bold mb-6">Ion Transport</h3>
        <div className="space-y-4">
          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">Ionic Conductivity</div>
            <div className="text-2xl font-quanto font-bold">
              {analysis.ionic_conductivity?.toExponential(2) || "N/A"} <span className="text-base">S/cm</span>
            </div>
          </div>

          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">Diffusion Coefficient</div>
            <div className="text-2xl font-quanto font-bold">
              {analysis.diffusion_coefficient?.toExponential(2) || "N/A"} <span className="text-base">cm²/s</span>
            </div>
          </div>

          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">Charge/Discharge Rate</div>
            <div className="text-2xl font-quanto font-bold text-purple-600 dark:text-purple-400">
              {analysis.charge_rate || "Fast"}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
