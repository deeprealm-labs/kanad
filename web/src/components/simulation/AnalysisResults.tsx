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

  // Enhanced molecular data
  const hasGeometry = results.geometry && Array.isArray(results.geometry) && results.geometry.length > 0;
  const hasOrbitalEnergies = results.orbital_energies && Array.isArray(results.orbital_energies) && results.orbital_energies.length > 0;
  const hasRDM1 = results.rdm1 && Array.isArray(results.rdm1) && results.rdm1.length > 0;
  const hasMolecularData = hasGeometry || hasOrbitalEnergies || hasRDM1 || results.nuclear_repulsion !== undefined;

  return (
    <div className="space-y-4 overflow-y-auto max-h-full">
      {/* Debug Info */}
      <div className="text-xs text-muted-foreground border border-dashed border-border rounded p-2">
        Analysis Sections: Energy({hasEnergyComponents ? "✓" : "✗"}), Bonding({hasBonding ? "✓" : "✗"}), Properties({hasProperties ? "✓" : "✗"}), Molecular Data({hasMolecularData ? "✓" : "✗"})
      </div>

      {/* Molecular Data (Enhanced) */}
      {hasMolecularData && (
        <div className="bg-card border border-border rounded-lg p-4">
          <h4 className="text-sm font-semibold mb-3 font-quando">
            Molecular Data
          </h4>
          <div className="space-y-3 text-xs">
            {/* Geometry */}
            {hasGeometry && (
              <div className="bg-muted rounded p-3">
                <div className="text-muted-foreground mb-2">
                  Molecular Geometry ({results.geometry.length} atoms)
                </div>
                <div className="font-mono space-y-1 max-h-32 overflow-y-auto">
                  {results.geometry.map((atom: [string, number[]], i: number) => (
                    <div key={i} className="flex justify-between items-center">
                      <span className="font-semibold text-blue-600 dark:text-blue-400 w-8">
                        {atom[0]}
                      </span>
                      <span className="text-muted-foreground flex-1 text-right">
                        ({atom[1][0].toFixed(4)}, {atom[1][1].toFixed(4)}, {atom[1][2].toFixed(4)})
                      </span>
                    </div>
                  ))}
                </div>
              </div>
            )}

            {/* Orbital Energies */}
            {hasOrbitalEnergies && (
              <div className="bg-muted rounded p-3">
                <div className="text-muted-foreground mb-2">
                  Orbital Energies ({results.orbital_energies.length} orbitals)
                </div>
                <div className="grid grid-cols-2 gap-2">
                  <div className="bg-green-50 dark:bg-green-900/20 border border-green-200 dark:border-green-800 rounded p-2">
                    <div className="text-green-700 dark:text-green-300 text-xs mb-1">HOMO</div>
                    <div className="font-mono font-semibold text-green-900 dark:text-green-100">
                      {results.orbital_energies[results.orbital_energies.length - 1].toFixed(4)} eV
                    </div>
                  </div>
                  {results.orbital_energies.length > 1 && (
                    <div className="bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-800 rounded p-2">
                      <div className="text-red-700 dark:text-red-300 text-xs mb-1">LUMO</div>
                      <div className="font-mono font-semibold text-red-900 dark:text-red-100">
                        {results.orbital_energies[0].toFixed(4)} eV
                      </div>
                    </div>
                  )}
                </div>
                {results.orbital_energies.length > 1 && (
                  <div className="mt-2 bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded p-2">
                    <div className="text-blue-700 dark:text-blue-300 text-xs mb-1">HOMO-LUMO Gap</div>
                    <div className="font-mono font-semibold text-blue-900 dark:text-blue-100">
                      {Math.abs(results.orbital_energies[0] - results.orbital_energies[results.orbital_energies.length - 1]).toFixed(4)} eV
                    </div>
                  </div>
                )}
              </div>
            )}

            {/* Density Matrix */}
            {hasRDM1 && (
              <div className="bg-muted rounded p-3">
                <div className="text-muted-foreground mb-2">
                  Reduced Density Matrix (RDM1)
                </div>
                <div className="font-mono font-semibold">
                  {results.rdm1.length} × {results.rdm1.length} matrix
                </div>
                <div className="text-xs text-muted-foreground mt-1">
                  Available for bonding and electronic structure analysis
                </div>
              </div>
            )}

            {/* Nuclear Repulsion */}
            {results.nuclear_repulsion !== undefined && (
              <div className="bg-muted rounded p-3">
                <div className="text-muted-foreground mb-1">Nuclear Repulsion Energy</div>
                <div className="font-mono font-semibold">
                  {results.nuclear_repulsion.toFixed(6)} Ha
                </div>
              </div>
            )}

            {/* Atom List */}
            {results.atoms && Array.isArray(results.atoms) && results.atoms.length > 0 && (
              <div className="bg-muted rounded p-3">
                <div className="text-muted-foreground mb-2">Atoms</div>
                <div className="flex flex-wrap gap-2">
                  {results.atoms.map((atom: string, i: number) => (
                    <span key={i} className="px-2 py-1 bg-blue-100 dark:bg-blue-900/30 text-blue-700 dark:text-blue-300 rounded font-semibold">
                      {atom}
                    </span>
                  ))}
                </div>
                <div className="text-xs text-muted-foreground mt-1">
                  {results.n_electrons} electrons, charge {results.charge || 0}, multiplicity {results.multiplicity || 1}
                </div>
              </div>
            )}
          </div>
        </div>
      )}

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

            {/* Individual Bond Analysis */}
            {bonding.bond_orders && bonding.bond_orders.bond_classification &&
              Object.keys(bonding.bond_orders.bond_classification).length > 0 && (
                <div className="bg-muted rounded p-3">
                  <div className="text-muted-foreground mb-2">Individual Bonds</div>
                  <div className="space-y-2">
                    {Object.entries(bonding.bond_orders.bond_classification).map(([bondKey, bondInfo]: [string, any]) => {
                      // Use bond_label if available, otherwise fall back to atom indices
                      const bondLabel = bondInfo.bond_label || `Atom ${bondInfo.atom_i} - Atom ${bondInfo.atom_j}`;

                      return (
                        <div key={bondKey} className="flex justify-between items-center border-b border-border pb-1 last:border-0">
                          <div className="font-medium">
                            {bondLabel}
                          </div>
                          <div className="flex gap-2 items-center">
                            <span className="text-blue-600 dark:text-blue-400 capitalize font-semibold">
                              {bondInfo.type}
                            </span>
                            <span className="font-mono text-muted-foreground text-xs">
                              ({bondInfo.order.toFixed(3)})
                            </span>
                          </div>
                        </div>
                      );
                    })}
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
                {!results.converged && (
                  <div className="text-xs text-amber-600 mt-1">
                    Try increasing max iterations in settings (recommended: 200+)
                  </div>
                )}
              </div>
            )}
          </div>
        </div>
      )}

      {/* Advanced Analysis Results */}
      {results.advanced_analysis && (
        <div className="bg-card border border-border rounded-lg p-4">
          <h4 className="text-sm font-semibold mb-3 font-quando">
            Advanced Analysis: {results.advanced_analysis.profile?.replace(/_/g, " ").replace(/\b\w/g, (l: string) => l.toUpperCase())}
          </h4>

          {results.advanced_analysis.status === "failed" ? (
            <div className="text-sm text-red-600 dark:text-red-400">
              {results.advanced_analysis.error || "Analysis failed"}
            </div>
          ) : (
            <div className="space-y-3">
              {Object.entries(results.advanced_analysis.results || {}).map(([moduleName, moduleResults]: [string, any]) => {
                const hasError = moduleResults?.error;
                const hasData = moduleResults && !hasError;

                return (
                  <div
                    key={moduleName}
                    className={`border rounded-lg p-3 ${
                      hasError
                        ? "border-red-300 bg-red-50 dark:border-red-800 dark:bg-red-900/20"
                        : "border-border bg-muted/50"
                    }`}
                  >
                    <h5 className="font-quando font-semibold text-sm mb-2 capitalize">
                      {moduleName.replace(/_/g, " ")}
                    </h5>

                    {hasError ? (
                      <p className="text-xs text-red-600 dark:text-red-400">{moduleResults.error}</p>
                    ) : moduleResults.message ? (
                      <div className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded p-3">
                        <p className="text-sm text-blue-700 dark:text-blue-300">
                          {moduleResults.message}
                        </p>
                        {moduleResults.future_features && Array.isArray(moduleResults.future_features) && (
                          <div className="mt-2">
                            <p className="text-xs text-blue-600 dark:text-blue-400 font-medium mb-1">Coming soon:</p>
                            <ul className="text-xs text-blue-600 dark:text-blue-400 list-disc list-inside space-y-1">
                              {moduleResults.future_features.map((feature: string, i: number) => (
                                <li key={i}>{feature}</li>
                              ))}
                            </ul>
                          </div>
                        )}
                      </div>
                    ) : (
                      <div className="space-y-2 text-xs">
                        {/* ADME Results */}
                        {moduleName === "adme" && moduleResults.adme_properties && (
                          <div className="grid grid-cols-2 gap-2">
                            <div className="bg-background rounded p-2">
                              <div className="text-muted-foreground mb-1">Molecular Weight</div>
                              <div className="font-semibold">
                                {moduleResults.descriptors?.molecular_weight?.toFixed(2) || "N/A"} g/mol
                              </div>
                            </div>
                            <div className="bg-background rounded p-2">
                              <div className="text-muted-foreground mb-1">LogP</div>
                              <div className="font-semibold">
                                {moduleResults.adme_properties?.logP?.toFixed(2) || "N/A"}
                              </div>
                            </div>
                            <div className="bg-background rounded p-2">
                              <div className="text-muted-foreground mb-1">BBB Permeability</div>
                              <div className="font-semibold">
                                {moduleResults.adme_properties?.bbb_permeability || "N/A"}
                              </div>
                            </div>
                            <div className="bg-background rounded p-2">
                              <div className="text-muted-foreground mb-1">Drug Likeness</div>
                              <div className="font-semibold">
                                {moduleResults.adme_properties?.drug_likeness || "N/A"}
                              </div>
                            </div>
                          </div>
                        )}

                        {/* Energy Analysis Results */}
                        {moduleName === "energy" && moduleResults.total_energy !== undefined && (
                          <div className="grid grid-cols-2 gap-2">
                            <div className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded p-2">
                              <div className="text-blue-700 dark:text-blue-300 mb-1">Total Energy</div>
                              <div className="font-mono font-semibold">
                                {moduleResults.total_energy.toFixed(6)} Ha
                              </div>
                            </div>
                            {moduleResults.electronic_energy !== undefined && (
                              <div className="bg-background rounded p-2">
                                <div className="text-muted-foreground mb-1">Electronic Energy</div>
                                <div className="font-mono font-semibold">
                                  {moduleResults.electronic_energy.toFixed(6)} Ha
                                </div>
                              </div>
                            )}
                          </div>
                        )}

                        {/* Frequencies Results */}
                        {moduleName === "frequencies" && moduleResults.frequencies && Array.isArray(moduleResults.frequencies) && (
                          <div className="space-y-2">
                            <div className="text-muted-foreground text-xs mb-2">
                              Vibrational Frequencies (cm⁻¹)
                            </div>
                            {moduleResults.frequencies.map((freq: number, i: number) => (
                              <div key={i} className="bg-background rounded p-2 flex justify-between items-center">
                                <span className="font-medium">Mode {i + 1}</span>
                                <span className="font-mono">{freq.toFixed(2)} cm⁻¹</span>
                              </div>
                            ))}
                            {moduleResults.zpe !== undefined && (
                              <div className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded p-2 mt-2">
                                <div className="text-blue-700 dark:text-blue-300 text-xs mb-1">Zero-Point Energy</div>
                                <div className="font-mono font-semibold">
                                  {moduleResults.zpe.toFixed(6)} Ha ({(moduleResults.zpe * 627.509).toFixed(2)} kcal/mol)
                                </div>
                              </div>
                            )}
                          </div>
                        )}

                        {/* Excited States Results */}
                        {moduleName === "excited_states" && moduleResults.excitation_energies && Array.isArray(moduleResults.excitation_energies) && (
                          <div className="space-y-2">
                            <div className="text-muted-foreground text-xs mb-2">
                              Electronic Excitations (Method: {moduleResults.method || "TDA"})
                            </div>
                            {moduleResults.excitation_energies.slice(0, 5).map((energy: number, i: number) => (
                              <div key={i} className="bg-background rounded p-2">
                                <div className="flex justify-between items-center mb-1">
                                  <span className="font-medium">S{i + 1}</span>
                                  <span className="font-mono">{energy.toFixed(4)} eV</span>
                                </div>
                                {moduleResults.wavelengths && moduleResults.wavelengths[i] && (
                                  <div className="flex justify-between items-center text-xs text-muted-foreground">
                                    <span>λ = {moduleResults.wavelengths[i].toFixed(2)} nm</span>
                                    {moduleResults.oscillator_strengths && moduleResults.oscillator_strengths[i] !== undefined && (
                                      <span>f = {moduleResults.oscillator_strengths[i].toFixed(4)}</span>
                                    )}
                                  </div>
                                )}
                              </div>
                            ))}
                          </div>
                        )}

                        {/* UV-Vis Results */}
                        {moduleName === "uv_vis" && moduleResults.excitation_energies && Array.isArray(moduleResults.excitation_energies) && (
                          <div className="space-y-2">
                            <div className="text-muted-foreground text-xs mb-2">
                              UV-Vis Absorption (Method: {moduleResults.method || "TDA"})
                            </div>
                            {moduleResults.excitation_energies.slice(0, 5).map((energy: number, i: number) => (
                              <div key={i} className="bg-background rounded p-2">
                                <div className="flex justify-between items-center mb-1">
                                  <span className="font-medium">S{i + 1}</span>
                                  <span className="font-mono">{energy.toFixed(4)} eV</span>
                                </div>
                                {moduleResults.wavelengths && moduleResults.wavelengths[i] && (
                                  <div className="flex justify-between items-center text-xs text-muted-foreground">
                                    <span>λ = {moduleResults.wavelengths[i].toFixed(2)} nm</span>
                                    {moduleResults.oscillator_strengths && moduleResults.oscillator_strengths[i] !== undefined && (
                                      <span>f = {moduleResults.oscillator_strengths[i].toFixed(4)}</span>
                                    )}
                                  </div>
                                )}
                              </div>
                            ))}
                          </div>
                        )}

                        {/* DOS (Density of States) Results */}
                        {moduleName === "dos" && moduleResults.energies && Array.isArray(moduleResults.energies) && (
                          <div className="space-y-2">
                            <div className="text-muted-foreground text-xs mb-2">
                              Density of States (Molecular Orbitals)
                            </div>
                            <div className="grid grid-cols-2 gap-2">
                              {moduleResults.homo_energy !== undefined && moduleResults.homo_energy !== null && (
                                <div className="bg-green-50 dark:bg-green-900/20 border border-green-200 dark:border-green-800 rounded p-2">
                                  <div className="text-green-700 dark:text-green-300 text-xs mb-1">HOMO</div>
                                  <div className="font-mono font-semibold text-green-900 dark:text-green-100">
                                    {moduleResults.homo_energy.toFixed(4)} eV
                                  </div>
                                </div>
                              )}
                              {moduleResults.lumo_energy !== undefined && moduleResults.lumo_energy !== null && (
                                <div className="bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-800 rounded p-2">
                                  <div className="text-red-700 dark:text-red-300 text-xs mb-1">LUMO</div>
                                  <div className="font-mono font-semibold text-red-900 dark:text-red-100">
                                    {moduleResults.lumo_energy.toFixed(4)} eV
                                  </div>
                                </div>
                              )}
                            </div>
                            {moduleResults.band_gap !== undefined && moduleResults.band_gap !== null && (
                              <div className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded p-2">
                                <div className="text-blue-700 dark:text-blue-300 text-xs mb-1">HOMO-LUMO Gap</div>
                                <div className="font-mono font-semibold text-blue-900 dark:text-blue-100">
                                  {moduleResults.band_gap.toFixed(4)} eV
                                </div>
                              </div>
                            )}
                            <div className="bg-background rounded p-2">
                              <div className="text-xs text-muted-foreground mb-1">
                                DOS Spectrum
                              </div>
                              <div className="text-xs">
                                {moduleResults.n_orbitals} molecular orbitals
                              </div>
                              <div className="text-xs mt-1">
                                Energy range: {Math.min(...moduleResults.energies).toFixed(2)} to {Math.max(...moduleResults.energies).toFixed(2)} eV
                              </div>
                              <div className="text-xs mt-1">
                                Broadening: {moduleResults.broadening} eV (Gaussian)
                              </div>
                            </div>
                          </div>
                        )}

                        {/* Thermochemistry Results */}
                        {moduleName === "thermochemistry" && moduleResults.temperature !== undefined && (
                          <div className="space-y-2">
                            <div className="text-muted-foreground text-xs mb-2">
                              Thermodynamic Properties at {moduleResults.temperature} K, {(moduleResults.pressure / 101325).toFixed(2)} atm
                            </div>
                            <div className="grid grid-cols-2 gap-2">
                              {moduleResults.e_elec !== undefined && (
                                <div className="bg-background rounded p-2">
                                  <div className="text-muted-foreground mb-1">Electronic Energy</div>
                                  <div className="font-mono font-semibold">
                                    {moduleResults.e_elec.toFixed(6)} Ha
                                  </div>
                                </div>
                              )}
                              {moduleResults.zpe !== undefined && (
                                <div className="bg-green-50 dark:bg-green-900/20 border border-green-200 dark:border-green-800 rounded p-2">
                                  <div className="text-green-700 dark:text-green-300 text-xs mb-1">Zero-Point Energy</div>
                                  <div className="font-mono font-semibold text-green-900 dark:text-green-100">
                                    {moduleResults.zpe.toFixed(6)} Ha
                                  </div>
                                </div>
                              )}
                              {moduleResults.h !== undefined && (
                                <div className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded p-2">
                                  <div className="text-blue-700 dark:text-blue-300 text-xs mb-1">Enthalpy (H)</div>
                                  <div className="font-mono font-semibold text-blue-900 dark:text-blue-100">
                                    {moduleResults.h.toFixed(6)} Ha
                                  </div>
                                </div>
                              )}
                              {moduleResults.s !== undefined && (
                                <div className="bg-purple-50 dark:bg-purple-900/20 border border-purple-200 dark:border-purple-800 rounded p-2">
                                  <div className="text-purple-700 dark:text-purple-300 text-xs mb-1">Entropy (S)</div>
                                  <div className="font-mono font-semibold text-purple-900 dark:text-purple-100">
                                    {moduleResults.s.toFixed(3)} cal/(mol·K)
                                  </div>
                                </div>
                              )}
                              {moduleResults.g !== undefined && (
                                <div className="bg-orange-50 dark:bg-orange-900/20 border border-orange-200 dark:border-orange-800 rounded p-2">
                                  <div className="text-orange-700 dark:text-orange-300 text-xs mb-1">Gibbs Free Energy (G)</div>
                                  <div className="font-mono font-semibold text-orange-900 dark:text-orange-100">
                                    {moduleResults.g.toFixed(6)} Ha
                                  </div>
                                </div>
                              )}
                              {moduleResults.cv !== undefined && (
                                <div className="bg-background rounded p-2">
                                  <div className="text-muted-foreground mb-1">Heat Capacity (Cv)</div>
                                  <div className="font-mono font-semibold">
                                    {moduleResults.cv.toFixed(3)} cal/(mol·K)
                                  </div>
                                </div>
                              )}
                              {moduleResults.cp !== undefined && (
                                <div className="bg-background rounded p-2">
                                  <div className="text-muted-foreground mb-1">Heat Capacity (Cp)</div>
                                  <div className="font-mono font-semibold">
                                    {moduleResults.cp.toFixed(3)} cal/(mol·K)
                                  </div>
                                </div>
                              )}
                              {moduleResults.e_thermal !== undefined && (
                                <div className="bg-background rounded p-2">
                                  <div className="text-muted-foreground mb-1">Thermal Correction</div>
                                  <div className="font-mono font-semibold">
                                    {moduleResults.e_thermal.toFixed(6)} Ha
                                  </div>
                                </div>
                              )}
                            </div>
                          </div>
                        )}

                        {/* Vibronic Results */}
                        {moduleName === "vibronic" && moduleResults.wavelengths && Array.isArray(moduleResults.wavelengths) && (
                          <div className="space-y-2">
                            <div className="text-muted-foreground text-xs mb-2">
                              Vibronic Spectrum (Method: {moduleResults.method || "Franck-Condon"})
                            </div>
                            {moduleResults.electronic_transition !== undefined && (
                              <div className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded p-2">
                                <div className="text-blue-700 dark:text-blue-300 text-xs mb-1">0-0 Transition</div>
                                <div className="font-mono font-semibold">
                                  {moduleResults.electronic_transition.toFixed(4)} eV ({(1239.84193 / moduleResults.electronic_transition).toFixed(2)} nm)
                                </div>
                              </div>
                            )}
                            <div className="bg-background rounded p-2">
                              <div className="text-xs text-muted-foreground mb-1">
                                Spectrum: {moduleResults.wavelengths.length} points
                              </div>
                              <div className="text-xs">
                                Range: {Math.min(...moduleResults.wavelengths).toFixed(0)} - {Math.max(...moduleResults.wavelengths).toFixed(0)} nm
                              </div>
                              {moduleResults.absorbance && Array.isArray(moduleResults.absorbance) && (
                                <div className="text-xs mt-1">
                                  Absorbance data computed ✓
                                </div>
                              )}
                              {moduleResults.emission && Array.isArray(moduleResults.emission) && (
                                <div className="text-xs mt-1">
                                  Emission data computed ✓
                                </div>
                              )}
                            </div>
                            {moduleResults.ground_frequencies && Array.isArray(moduleResults.ground_frequencies) && (
                              <div className="bg-background rounded p-2">
                                <div className="text-xs text-muted-foreground mb-1">
                                  Ground State Frequencies
                                </div>
                                <div className="text-xs font-mono">
                                  {moduleResults.ground_frequencies.slice(0, 3).map((f: number) => f.toFixed(1)).join(', ')} cm⁻¹
                                  {moduleResults.ground_frequencies.length > 3 && ` (+${moduleResults.ground_frequencies.length - 3} more)`}
                                </div>
                              </div>
                            )}
                          </div>
                        )}
                      </div>
                    )}
                  </div>
                );
              })}
            </div>
          )}
        </div>
      )}
    </div>
  );
}
