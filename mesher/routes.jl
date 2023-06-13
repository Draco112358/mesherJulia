using Genie.Router, Genie.Renderer, Genie.Renderer.Html, Genie.Renderer.Json, Genie.Requests, JSON

include("lib/mesher.jl")

route("/") do
  serve_static_file("welcome.html")
end


route("/meshing" ,method="POST") do 
    return JSON.json(doMeshing(jsonpayload()))
end