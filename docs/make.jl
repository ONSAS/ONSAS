
using Documenter, ONSAS_docs

DocMeta.setdocmeta!(ONSAS_docs, :DocTestSetup,
  :(using ONSAS_docs); recursive=true)

makedocs(
  sitename="ONSAS",
  format=Documenter.HTML(prettyurls=get(ENV, "CI", nothing) == "true"),
  pages=[
    "Home" => "index.md",
    "Learning by examples" => Any["Static Von-Mises Truss"=>"examples/staticVonMisesTruss.md",
      "Spring-mass system"=>"examples/springMass.md",
      "Uniaxial extension"=>"examples/uniaxialExtension.md",
      "Uniaxial compression"=>"examples/uniaxialCompression.md",
      "Cantilever Beam"=>"examples/cantileverBeam.md",
      "Linear beam vibration"=>"examples/beamLinearVibration.md",
      "Plane strain ring"=>"examples/ringPlaneStrain.md",
      # "Linear aerodynamics"     => "examples/linearAerodynamics.md",
      "Reconfiguration beam"=>"examples/dragBeamReconfiguration.md",
      "Propeller model"=>"examples/simplePropeller.md"],
      # "Non-linear aerodynamics"=>"examples/nonLinearAerodynamics.md"],
      "User guide" => Any["Installation"=>"install.md",
      "Creating Models"=>"creatingModels.md",
      "References"=>"references.md"]
    #"Developer guide" => Any["ONSAS_solve" =>  "ONSAS_solve.md"]
  ]
)

deploydocs(
  repo="github.com/ONSAS/ONSAS.git",
  push_preview=true)
