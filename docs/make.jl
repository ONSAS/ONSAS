
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
      "Learning by examples" => Any["Static Von-Mises Truss" =>  "staticVonMisesTruss.md",
                                    "Uniaxial extension" => "uniaxialExtension.md",
                                    "Cantilever Beam" => "cantileverBeam.md",
                                    "Solid with inclusion" => "semiSphereWithInclusion.md"],
      "User guide" => Any["Installation" =>  "howtouse/install.md",
                          "Creating Models"  => "howtouse/creatingModels.md",
                          "Corotational frame"  => "corotationalFrameElement.md",
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
