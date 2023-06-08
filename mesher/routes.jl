using Genie.Router

route("/") do
  serve_static_file("welcome.html")
end


route("/meshing" ,method="POST") do 
    return JSON.json(doMeshing(jsonpayload()["mesherInput"]))
end