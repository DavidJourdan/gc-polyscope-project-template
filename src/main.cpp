#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/stripe_patterns.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"

#include "args/args.hxx"
#include "imgui.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

// Some algorithm parameters
float param1 = 42.0;

// Example computation function -- this one computes and registers a scalar
// quantity
void doWork() {
  polyscope::warning("Computing Gaussian curvature.\nalso, parameter value = " +
                     std::to_string(param1));

  geometry->requireVertexGaussianCurvatures();
  psMesh->addVertexScalarQuantity("curvature",
                                  geometry->vertexGaussianCurvatures,
                                  polyscope::DataType::SYMMETRIC);
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {

  if (ImGui::Button("do work")) {
    doWork();
  }

  ImGui::SliderFloat("param", &param1, 0., 100.);
}

int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("geometry-central & Polyscope example project");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help &h) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load mesh
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(args::get(inputFilename));

  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath(args::get(inputFilename)),
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  // Set vertex tangent spaces
  geometry->requireVertexTangentBasis();
  VertexData<Vector3> vBasisX(*mesh);
  for (Vertex v : mesh->vertices()) {
    vBasisX[v] = geometry->vertexTangentBasis[v][0];
  }
  psMesh->setVertexTangentBasisX(vBasisX);

  auto vField =
      geometrycentral::surface::computeCurvatureAlignedVertexDirectionField(*geometry, 2);
  psMesh->addVertexIntrinsicVectorQuantity("VF", vField, 2);

  geometry->requireEdgeLengths();
  double avgLength = geometry->edgeLengths.toVector().sum() / mesh->nEdges();

  VertexData<double> frequencies(*mesh, 2 / avgLength);
  const auto& [stripeValues, stripeIndices, fieldIndices] = computeStripePattern(*geometry, frequencies, vField);
  const auto& [vertices, edges] = extractPolylinesFromStripePattern(*geometry, stripeValues, stripeIndices, fieldIndices, vField, true);

  std::vector<std::pair<size_t, int>> stripeCount;
  for (size_t iF = 0; iF < mesh->nFaces(); iF++) {
    if (stripeIndices[iF] != 0) {
      stripeCount.push_back(std::make_pair(iF, stripeIndices[iF]));
    } 
  }
  psMesh->addFaceCountQuantity("Stripe indices", stripeCount);

  std::vector<std::pair<size_t, int>> fieldCount;
  for (size_t iF = 0; iF < mesh->nFaces(); iF++) {
    if (fieldIndices[iF] != 0) {
      fieldCount.push_back(std::make_pair(iF, fieldIndices[iF]));
    } 
  }
  psMesh->addFaceCountQuantity("Field indices", fieldCount);

  polyscope::registerCurveNetwork("stripes", vertices, edges);

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
