from flask import Flask, jsonify, make_response, request
from flask_cors import cross_origin
from mesher_py.mesher import doMeshing
from solver_py.solver import doSolving

app = Flask(__name__)


@app.route("/")
def hello_from_root():
    return jsonify(message='Hello from root!')


@app.route("/hello")
def hello():
    return jsonify(message='Hello from path!')

@app.route("/meshing", methods=['POST'])
@cross_origin()
def meshing():
    inputMesher = request.get_json()
    return doMeshing(inputMesher)

@app.route("/solving", methods=['POST'])
@cross_origin()
def solving():
    mesherOutput = request.get_json()['mesherOutput']
    solverInput = request.get_json()['solverInput']
    return doSolving(mesherOutput, solverInput)


@app.errorhandler(404)
def resource_not_found(e):
    return make_response(jsonify(error='Not found!'), 404)
