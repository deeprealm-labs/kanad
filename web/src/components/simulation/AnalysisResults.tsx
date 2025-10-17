"use client";

interface AnalysisResultsProps {
  results: any;
}

export default function AnalysisResults({ results }: AnalysisResultsProps) {
  // Debug logging
  console.log("AnalysisResults received results:", results);
  console.log("Has analysis?", !!results?.analysis);
  if (results?.analysis) {
    console.log("Analysis keys:", Object.keys(results.analysis));
    console.log("Energy components:", results.analysis.energy_components);
    console.log("Bonding:", results.analysis.bonding);
    console.log("Properties:", results.analysis.properties);
  }

  if (!results || !results.analysis) {
    return (
      <div className="text-sm text-muted-foreground text-center py-8">
        No analysis data available
      </div>
    );
  }

  const { analysis } = results;
  const energyComponents = analysis.energy_components || {};
  const bonding = analysis.bonding || {};
  const properties = analysis.properties || {};

  const hasEnergyComponents = Object.keys(energyComponents).length > 0;
  const hasBonding = bonding.bond_type !== undefined;
  const hasProperties = Object.keys(properties).length > 0;

  return (
    <div className="space-y-4 overflow-y-auto max-h-full">
      {/* Debug Info */}
      <div className="text-xs text-muted-foreground border border-dashed border-border rounded p-2">
        Analysis Sections: Energy({hasEnergyComponents ? "✓" : "✗"}), Bonding({hasBonding ? "✓" : "✗"}), Properties({hasProperties ? "✓" : "✗"})
      </div>

      {/* Energy Components */}
      {hasEnergyComponents && (
        <div className="bg-card border border-border rounded-lg p-4">
          <h4 className="text-sm font-semibold mb-3 font-quando">
            Energy Components
          </h4>
          <div className="grid grid-cols-2 gap-3 text-xs">
            {energyComponents.nuclear_repulsion !== undefined && (
              <div className="bg-muted rounded p-2">
                <div className="text-muted-foreground mb-1">Nuclear Repulsion</div>
                <div className="font-mono font-semibold">
                  {energyComponents.nuclear_repulsion.toFixed(6)} Ha
                </div>
              </div>
            )}
            {energyComponents.one_electron !== undefined && (
              <div className="bg-muted rounded p-2">
                <div className="text-muted-foreground mb-1">One-Electron</div>
                <div className="font-mono font-semibold">
                  {energyComponents.one_electron.toFixed(6)} Ha
                </div>
              </div>
            )}
            {energyComponents.coulomb !== undefined && (
              <div className="bg-muted rounded p-2">
                <div className="text-muted-foreground mb-1">Coulomb</div>
                <div className="font-mono font-semibold">
                  {energyComponents.coulomb.toFixed(6)} Ha
                </div>
              </div>
            )}
            {energyComponents.exchange !== undefined && (
              <div className="bg-muted rounded p-2">
                <div className="text-muted-foreground mb-1">Exchange</div>
                <div className="font-mono font-semibold">
                  {energyComponents.exchange.toFixed(6)} Ha
                </div>
              </div>
            )}
            {energyComponents.two_electron !== undefined && (
              <div className="bg-muted rounded p-2">
                <div className="text-muted-foreground mb-1">Two-Electron</div>
                <div className="font-mono font-semibold">
                  {energyComponents.two_electron.toFixed(6)} Ha
                </div>
              </div>
            )}
            {energyComponents.total !== undefined && (
              <div className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded p-2">
                <div className="text-blue-700 dark:text-blue-300 mb-1">Total Energy</div>
                <div className="font-mono font-bold text-blue-900 dark:text-blue-100">
                  {energyComponents.total.toFixed(6)} Ha
                </div>
              </div>
            )}
          </div>
        </div>
      )}

      {/* Bonding Analysis */}
      {bonding.bond_type && (
        <div className="bg-card border border-border rounded-lg p-4">
          <h4 className="text-sm font-semibold mb-3 font-quando">
            Bonding Analysis
          </h4>
          <div className="space-y-3 text-xs">
            <div className="bg-muted rounded p-3">
              <div className="text-muted-foreground mb-1">Bond Type</div>
              <div className="font-semibold capitalize">
                {bonding.bond_type.bonding_type || "N/A"}
              </div>
            </div>

            {bonding.bond_type.homo_lumo_gap !== undefined && (
              <div className="grid grid-cols-2 gap-3">
                <div className="bg-muted rounded p-2">
                  <div className="text-muted-foreground mb-1">HOMO-LUMO Gap</div>
                  <div className="font-mono font-semibold">
                    {bonding.bond_type.homo_lumo_gap.toFixed(6)} Ha
                  </div>
                </div>
                {bonding.bond_type.homo_lumo_gap_ev !== undefined && (
                  <div className="bg-muted rounded p-2">
                    <div className="text-muted-foreground mb-1">HOMO-LUMO Gap</div>
                    <div className="font-mono font-semibold">
                      {bonding.bond_type.homo_lumo_gap_ev.toFixed(4)} eV
                    </div>
                  </div>
                )}
              </div>
            )}

            {bonding.bond_type.characteristics &&
              Array.isArray(bonding.bond_type.characteristics) &&
              bonding.bond_type.characteristics.length > 0 && (
                <div className="bg-muted rounded p-3">
                  <div className="text-muted-foreground mb-2">Characteristics</div>
                  <ul className="list-disc list-inside space-y-1">
                    {bonding.bond_type.characteristics.map((char: string, i: number) => (
                      <li key={i} className="text-foreground">
                        {char}
                      </li>
                    ))}
                  </ul>
                </div>
              )}

            {bonding.bond_orders && bonding.bond_orders.bond_orders && (
              <div className="bg-muted rounded p-3">
                <div className="text-muted-foreground mb-2">Bond Orders</div>
                <div className="font-mono text-xs">
                  {JSON.stringify(bonding.bond_orders.bond_orders, null, 2)}
                </div>
              </div>
            )}
          </div>
        </div>
      )}

      {/* Molecular Properties */}
      {Object.keys(properties).length > 0 && (
        <div className="bg-card border border-border rounded-lg p-4">
          <h4 className="text-sm font-semibold mb-3 font-quando">
            Molecular Properties
          </h4>
          <div className="space-y-3 text-xs">
            {properties.dipole_moment !== undefined && (
              <div className="bg-muted rounded p-3">
                <div className="text-muted-foreground mb-1">Dipole Moment</div>
                <div className="font-mono font-semibold text-base">
                  {properties.dipole_moment.toFixed(6)} D
                </div>
              </div>
            )}

            {properties.dipole_vector && Array.isArray(properties.dipole_vector) && (
              <div className="bg-muted rounded p-3">
                <div className="text-muted-foreground mb-2">Dipole Vector</div>
                <div className="font-mono space-y-1">
                  <div>x: {properties.dipole_vector[0]?.toFixed(6) || "0.000000"}</div>
                  <div>y: {properties.dipole_vector[1]?.toFixed(6) || "0.000000"}</div>
                  <div>z: {properties.dipole_vector[2]?.toFixed(6) || "0.000000"}</div>
                </div>
              </div>
            )}

            {properties.center_of_mass && Array.isArray(properties.center_of_mass) && (
              <div className="bg-muted rounded p-3">
                <div className="text-muted-foreground mb-2">Center of Mass</div>
                <div className="font-mono space-y-1">
                  <div>x: {properties.center_of_mass[0]?.toFixed(6) || "0.000000"}</div>
                  <div>y: {properties.center_of_mass[1]?.toFixed(6) || "0.000000"}</div>
                  <div>z: {properties.center_of_mass[2]?.toFixed(6) || "0.000000"}</div>
                </div>
              </div>
            )}

            {properties.center_of_charge && Array.isArray(properties.center_of_charge) && (
              <div className="bg-muted rounded p-3">
                <div className="text-muted-foreground mb-2">Center of Charge</div>
                <div className="font-mono space-y-1">
                  <div>x: {properties.center_of_charge[0]?.toFixed(6) || "0.000000"}</div>
                  <div>y: {properties.center_of_charge[1]?.toFixed(6) || "0.000000"}</div>
                  <div>z: {properties.center_of_charge[2]?.toFixed(6) || "0.000000"}</div>
                </div>
              </div>
            )}
          </div>
        </div>
      )}

      {/* Additional Results */}
      {results.hf_energy !== undefined && (
        <div className="bg-card border border-border rounded-lg p-4">
          <h4 className="text-sm font-semibold mb-3 font-quando">
            Additional Results
          </h4>
          <div className="grid grid-cols-2 gap-3 text-xs">
            <div className="bg-muted rounded p-2">
              <div className="text-muted-foreground mb-1">HF Energy</div>
              <div className="font-mono font-semibold">
                {results.hf_energy.toFixed(6)} Ha
              </div>
            </div>
            {results.correlation_energy !== undefined && (
              <div className="bg-muted rounded p-2">
                <div className="text-muted-foreground mb-1">Correlation Energy</div>
                <div className="font-mono font-semibold">
                  {results.correlation_energy.toFixed(6)} Ha
                </div>
              </div>
            )}
            {results.iterations !== undefined && (
              <div className="bg-muted rounded p-2">
                <div className="text-muted-foreground mb-1">Iterations</div>
                <div className="font-mono font-semibold">
                  {results.iterations}
                </div>
              </div>
            )}
            {results.converged !== undefined && (
              <div className="bg-muted rounded p-2">
                <div className="text-muted-foreground mb-1">Converged</div>
                <div className={`font-semibold ${results.converged ? "text-green-600" : "text-red-600"}`}>
                  {results.converged ? "Yes" : "No"}
                </div>
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  );
}
