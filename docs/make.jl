
using Documenter, ONSAS_docs

DocMeta.setdocmeta!(ONSAS_docs, :DocTestSetup,
                    :(using ONSAS_docs); recursive=true)

makedocs(
  sitename = "ONSAS.m",
  format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
  pages = [
    "Home" => "index.md",
    "Learning by examples" => Any["Static Von-Mises Truss"  =>  "examples/staticVonMisesTruss.md",
                                  "Spring-mass system"      => "examples/springMass.md",
                                  "Uniaxial extension"      => "examples/uniaxialExtension.md",
                                  "Cantilever Beam"         => "examples/cantileverBeam.md",
                                  "Linear beam vibration"   => "examples/beamLinearVibration.md",
                                  "Linear elastic cylinder" => "examples/linearCylinderPlaneStrain.md",
                                  "Linear aerodynamics"     => "examples/linearAerodynamics.md",
                                  "Non-linear aerodynamics" => "examples/nonLinearAerodynamics.md"],
    "User guide" => Any["Installation" =>  "howtouse/install.md",
                        "Creating Models"  => "howtouse/creatingModels.md",
                        "References"  => "theory/references.md"]
    #"Developer guide" => Any["ONSAS_solve" =>  "ONSAS_solve.md"]
  ],
  strict = false
)
  
deploydocs(
  repo = "github.com/ONSAS/ONSAS.m.git",
  push_preview=true )
