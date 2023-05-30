using Genie.Router

route("/") do
  serve_static_file("welcome.html")
end

route("/hello") do
  return hello()
end

route("/meshing" ,method="POST") do 
    # return JSON.json(doMeshing(jsonpayload()["mesherOutput"], jsonpayload()["solverInput"], jsonpayload()["solverAlgoParams"], client))
    return "meshing"
end