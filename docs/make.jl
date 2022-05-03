
using Documenter, ONSAS_docs

DocMeta.setdocmeta!(ONSAS_docs, :DocTestSetup,
                    :(using ONSAS_docs); recursive=true)

# sets default format or read provided argument
if length(ARGS)==0
  outputFormat = "html"
else
  outputFormat = ARGS[1]
end

# sets different format and pages settings for each output format
if outputFormat == "pdf"
  makedocs(
    sitename = "ONSAS.m",
    modules = [ONSAS_docs],
    format = Documenter.LaTeX(),
    pages = [
      "index.md",
      "uniaxialExtension.md"
    ]
  )

elseif outputFormat == "html"

  makedocs(
    sitename = "ONSAS.m",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
      "Home" => "index.md",
      "Learning by examples" => Any["Static Von-Mises Truss" =>  "examples/staticVonMisesTruss.md",
                                     "Spring-mass system" => "examples/springMass.md",
                                    "Uniaxial extension" => "examples/uniaxialExtension.md",
                                    "Cantilever Beam" => "examples/cantileverBeam.md",
                                    "Solid with inclusion" => "examples/semiSphereWithInclusion.md",
                                    "Linear beam vibration" => "examples/beamLinearVibration.md",
                                    "Linear aerodynamics" => "examples/linearAerodynamics.md",
                                    "Non-linear aerodynamics" => "examples/nonLinearAerodynamics.md"],
      "User guide" => Any["Installation" =>  "howtouse/install.md",
                          "Creating Models"  => "howtouse/creatingModels.md",
#                          "Corotational frame"  => "corotationalFrameElement.md",
                          "References"  => "theory/references.md"]
      #"Developer guide" => Any["ONSAS_solve" =>  "ONSAS_solve.md"]
    ],
    strict = false
  )

  deploydocs(
      repo = "github.com/ONSAS/ONSAS.m.git",
      push_preview=true
  )

end
#

#         "Guide" => Any["Installation"     => "howtouse/install.md",
#         "Theory" => Any["Virtual mechanical work " => "theory/prinMechWork.md",
#                         "References"               => "theory/references.md"  ],
#         "About" => "about.md",
#     ],
#     strict = false
# )
#
