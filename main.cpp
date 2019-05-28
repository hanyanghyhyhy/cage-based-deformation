// DF: Deformation/ Deformation Cage
// CF: Control Point/ Control Vertex
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/edges.h>
#include <igl/file_exists.h>

// Threshold to avoid numerical instability
double EPSILON = 1.0e-19;
//Eigen::MatrixXd headPointsCage, tailPointsCage, headPointsCageDF, tailPointsCageDF;
Eigen::MatrixXd mvcWeights1,mvcWeights2,mvcWeights3;

// declare the function to assign weight
Eigen::MatrixXd operator*(const Eigen::MatrixXd & weights, Eigen::MatrixXd & cageVertex);

// function to compute determinant for matrix U which is composed of three direction unit vectors
inline double computeDeterminant(Eigen::MatrixXd u0, Eigen::MatrixXd u1, Eigen::MatrixXd u2) {
    double determinant = u0(0) * (u1(1) * u2(2) - u2(1) * u1(2))
                         - u0(1) * (u1(0) * u2(2) - u2(0) * u1(2))
                         + u0(2) * (u1(0) * u2(1) - u2(0) * u1(1));
    return determinant;
}

class TargetModel {
public:
    // constructor
    TargetModel(std::string fileName) : nv_len(0.5), point_size(8), line_width(0.5), sel_vidx(7), mode(0) {
        //initial vertices and faces
        if (!igl::file_exists(fileName)) {
            std::cout << "[error] cannot locate model file at " << fileName << "\nPress any key to exit\n";
            char c;
            std::cin >> c;
            exit(1);
        }

        // read mesh model
        igl::readOFF(fileName, m_V, m_F);
        // compute face normals
        igl::per_face_normals(m_V, m_F, m_FN);
    }
    // deconstructor
    ~TargetModel() {}

    // model vertices, faces and face normals
    Eigen::MatrixXd m_V;
    Eigen::MatrixXi m_F;
    Eigen::MatrixXd m_FN;

    // parameters for UI settings
    float nv_len;
    float point_size;
    float line_width;
    int sel_vidx;
    int mode;


    void implementNoise(float noiseLevel) {
        // compute the bounding box dimensions
        Eigen::Vector3d minBound = m_V.colwise().minCoeff();
        Eigen::Vector3d maxBound = m_V.colwise().maxCoeff();
        // add noise to each column
        Eigen::MatrixXd noise = Eigen::MatrixXd::Zero(m_V.rows(), m_V.cols());
        for (int i = 0; i < m_V.cols(); i++) {
            noise.col(i) = noiseLevel * (maxBound(i) - minBound(i)) * (Eigen::MatrixXd::Random(m_V.rows(), 1));
        }
        m_V += noise;
    }


private:
};

class ControlCage {
public:
    // constructor
    ControlCage(std::string fileName) {
        if (!igl::file_exists(fileName)) {
            std::cout << "[error] cannot locate model file at " << fileName << "\nPress any key to exit\n";
            char c;
            std::cin >> c;
            exit(1);
        }
        igl::readOFF(fileName, m_V, m_F);
    }
    // deconstructor
    ~ControlCage() {}

    Eigen::MatrixXd m_V, headPoints, tailPoints;
    Eigen::MatrixXi m_F;

    void loadCage() {
        // The edges function returns the head and the trail position of the edges
        // use the position vectors to retrieve the edges of the cage
        // store them in headPoints and tailPoints.
        Eigen::MatrixXi Edges;
        igl::edges(m_F, Edges);
        headPoints = Eigen::MatrixXd(Edges.rows(), 3);
        tailPoints = Eigen::MatrixXd(Edges.rows(), 3);

        for (int i = 0; i < Edges.rows(); i++) {
            headPoints.row(i) = m_V.row(Edges(i, 0));
            tailPoints.row(i) = m_V.row(Edges(i, 1));
        }
    }

private:
};




// Main function starts
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {

    std::cout << "Key: " << key << " " << (unsigned int) key << std::endl;
    if (key == 'q' || key == 'Q') {
        exit(0);
    }
    return false;
}

// x lies on t, use 2D barycentric coordinates
Eigen::MatrixXd computeBarycentricCoordinates(int index, int rows, Eigen::MatrixXi controlFace, double theta0, double theta1,
                              double theta2, double l0, double l1, double l2) {

    Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(rows, 1);
    // use sub-triangle area to represent weights
    weights(controlFace(index, 0)) = sin(theta0) * l1 * l2;
    weights(controlFace(index, 1)) = sin(theta1) * l0 * l2;
    weights(controlFace(index, 2)) = sin(theta2) * l0 * l1;

    return weights;
}

// This function only compute the weight between single model point and control point.
Eigen::MatrixXd computeSingleWeight(Eigen::MatrixXd modelVertex, Eigen::MatrixXd controlVertex, Eigen::MatrixXi controlFace) {

    // initialise the weight vector corresponding to current model vertex (1, Number of cage vertex)
    Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(1, controlVertex.rows());

    // checks the distance between control point and model point
    for (int i = 0; i < controlVertex.rows(); i++) {
        // extract a single p_j
        Eigen::MatrixXd individualCP = controlVertex.row(i);
        // d_i = |p_j - x|
        Eigen::Vector3d d(modelVertex(0) - individualCP(0), modelVertex(1) - individualCP(1),
                          modelVertex(2) - individualCP(2));

        // if the distance between control vertex and model vertex is too close, return 1 to weight.
        // fully weight this vertex and others are zero
        if (d.norm() <= EPSILON) {
            weights(i) = 1.0;
            return weights;
        }
    }

    for (int j = 0; j < controlFace.rows(); j++) {

        // for each triangle face, extract its vertices
        Eigen::MatrixXd p0 = controlVertex.row(controlFace(j, 0));
        Eigen::MatrixXd p1 = controlVertex.row(controlFace(j, 1));
        Eigen::MatrixXd p2 = controlVertex.row(controlFace(j, 2));

        // compute u_i = (p_j - x) / ||p_j - x||  Normalised direction vector
        Eigen::Vector3d d0(modelVertex(0) - p0(0), modelVertex(1) - p0(1), modelVertex(2) - p0(2));
        Eigen::Vector3d d1(modelVertex(0) - p1(0), modelVertex(1) - p1(1), modelVertex(2) - p1(2));
        Eigen::Vector3d d2(modelVertex(0) - p2(0), modelVertex(1) - p2(1), modelVertex(2) - p2(2));
        Eigen::Vector3d u0 = d0.normalized();
        Eigen::Vector3d u1 = d1.normalized();
        Eigen::Vector3d u2 = d2.normalized();

        // For the spherical triangle: compute the edge length of the triangle which has the same vertices with the spherical triangle
        double l0 = (u1 - u2).norm();
        double l1 = (u0 - u2).norm();
        double l2 = (u0 - u1).norm();
        // arc length of the spherical triangle
        double theta0 = 2 * asin(l0 * 0.5);
        double theta1 = 2 * asin(l1 * 0.5);
        double theta2 = 2 * asin(l2 * 0.5);
        // half of the perimeter of the spherical triangle
        double h = (theta0 + theta1 + theta2) * 0.5;
        // check whether the model vertex lies on the plane containing the spherical triangle
        if ((M_PI - h) <= EPSILON) {
            // if so, use 2D barycentric coordinates to compute weight
            return computeBarycentricCoordinates(j, controlVertex.rows(), controlFace, theta0, theta1, theta2, l0, l1,
                                                 l2);
        }

        // compute the cosine of dihedral angles between side faces
        double c0 = (2 * sin(h) * sin(h - theta0)) / (sin(theta1) * sin(theta2)) - 1.0;
        double c1 = (2 * sin(h) * sin(h - theta1)) / (sin(theta0) * sin(theta2)) - 1.0;
        double c2 = (2 * sin(h) * sin(h - theta2)) / (sin(theta0) * sin(theta1)) - 1.0;
        // check the range of the angle
        double determinant = computeDeterminant(u0, u1, u2);
        double sign = (determinant < EPSILON) ? -1 : 1;
        // compute the sine of dihedral angles between side faces
        double s0 = sign * determinant * sqrt(1.0 - c0 * c0);
        double s1 = sign * determinant * sqrt(1.0 - c1 * c1);
        double s2 = sign * determinant * sqrt(1.0 - c2 * c2);

        // x lies outside t on the same plane, ignore t
        if (s0 <= EPSILON || s1 <= EPSILON || s2 <= EPSILON) {
            continue;
        }
        // compute weight and add them to the sum
        weights(controlFace(j, 0)) += (theta0 - c1 * theta2 - c2 * theta1) / (2 * d0.norm() * sin(theta2) * sqrt(1 - c1 * c1));
        weights(controlFace(j, 1)) += (theta1 - c0 * theta2 - c2 * theta0) / (2 * d1.norm() * sin(theta0) * sqrt(1 - c2 * c2));
        weights(controlFace(j, 2)) += (theta2 - c1 * theta0 - c0 * theta1) / (2 * d2.norm() * sin(theta1) * sqrt(1 - c0 * c0));

    }

    return weights;
}

// Compute the weights between model vertices and control cage before deformation
// This function is called before the deformation to expedite the computation
// For each vertex on the model surface, we compute the weights between it and all the control points
Eigen::MatrixXd computeWeight(Eigen::MatrixXd modelVertex, Eigen::MatrixXd cageVertex, Eigen::MatrixXi cageFace) {
    // the size of weight matrix is (Number of model vertex, Number of cage vertex)
    // for each model vertex, its weight is evaluated on all cage vertices
    Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(modelVertex.rows(), cageVertex.rows());
    // go through each vertex in the model and compute the corresponding weight
    for (int i = 0; i < modelVertex.rows(); i++) {
        weights.row(i) = computeSingleWeight(modelVertex.row(i), cageVertex, cageFace);
    }

    return weights;
}

// After deformation, update the vetices in model using computed weights.
//Eigen::MatrixXd updateDeformedVertex(Eigen::MatrixXd weights, Eigen::MatrixXd cageVertex) {
//    Eigen::MatrixXd updatedVertex = Eigen::MatrixXd::Zero(weights.rows(), 3);
//
//    for (int i = 0; i < weights.rows(); i++) {
//        Eigen::MatrixXd singleRow = Eigen::MatrixXd::Zero(1, 3);
//        double weightSum = 0.0;
//        for (int j = 0; j < weights.cols(); j++) {
//            singleRow += weights(i, j) * cageVertex.row(j);
//            weightSum += weights(i, j);
//        }
//        updatedVertex.row(i) = singleRow / weightSum;
//    }
//
//    return updatedVertex;
//}

// Use operator overload
// assign the vertex weight on the deformed cage to get the deformed model
Eigen::MatrixXd operator*(const Eigen::MatrixXd &weights, Eigen::MatrixXd &cageVertex) {
    // initialise the deformed model
    Eigen::MatrixXd updatedVertex = Eigen::MatrixXd::Zero(weights.rows(), 3);
    // assign weight
    for (int i = 0; i < weights.rows(); i++) {
        Eigen::MatrixXd singleRow = Eigen::MatrixXd::Zero(1, 3);
        double weightSum = 0.0;
        for (int j = 0; j < weights.cols(); j++) {
            singleRow += weights(i, j) * cageVertex.row(j);
            weightSum += weights(i, j);
        }
        // normalise the weight
        updatedVertex.row(i) = singleRow / weightSum;
    }

    return updatedVertex;
}

int main(int argc, char *argv[]) {
    // model path
    std::string modelPath1 = "../meshes/beast.off";
    std::string cagePath1 = "../meshes/beast_cage.off";
    std::string dfPath1 = "../meshes/beast_cage_deformed.off";
    std::string modelPath2 = "../meshes/bench.off";
    std::string cagePath2 = "../meshes/bench_cage.off";
    std::string dfPath2 = "../meshes/bench_cage_deformed.off";
    std::string modelPath3 = "../meshes/cactus.off";
    std::string cagePath3 = "../meshes/cactus_cage.off";
    std::string dfPath3 = "../meshes/cactus_cage_deformed.off";

    // load model 1
    TargetModel targetModel1 = TargetModel(modelPath1);
    ControlCage controlCage1 = ControlCage(cagePath1);
    // Calculate the weight in advance to exoedite the computation
    mvcWeights1 = computeWeight(targetModel1.m_V, controlCage1.m_V, controlCage1.m_F);
    ControlCage deformedCage1 = ControlCage(dfPath1);
    Eigen::MatrixXd updatedVertex1 = mvcWeights1 * deformedCage1.m_V;
    deformedCage1.loadCage();
    controlCage1.loadCage();

    // load model 2
    TargetModel targetModel2 = TargetModel(modelPath2);
    ControlCage controlCage2 = ControlCage(cagePath2);
    // Calculate the weight in advance to exoedite the computation
    mvcWeights2 = computeWeight(targetModel2.m_V, controlCage2.m_V, controlCage2.m_F);
    ControlCage deformedCage2 = ControlCage(dfPath2);
    Eigen::MatrixXd updatedVertex2 = mvcWeights2 * deformedCage2.m_V;
    deformedCage2.loadCage();
    controlCage2.loadCage();

    // load model 3
    TargetModel targetModel3 = TargetModel(modelPath3);
    ControlCage controlCage3 = ControlCage(cagePath3);
    // Calculate the weight in advance to exoedite the computation
    mvcWeights3 = computeWeight(targetModel3.m_V, controlCage3.m_V, controlCage3.m_F);
    ControlCage deformedCage3 = ControlCage(dfPath3);
    Eigen::MatrixXd updatedVertex3= mvcWeights3 * deformedCage3.m_V;
    deformedCage3.loadCage();
    controlCage3.loadCage();


    // Init the viewer
    igl::opengl::glfw::Viewer viewer;
    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    menu.callback_draw_viewer_menu = [&]() {
        // Add new group
        if (ImGui::CollapsingHeader("Deformation Model 1", ImGuiTreeNodeFlags_DefaultOpen)) {
            // Add buttons

            if (ImGui::Button("Show / Reset Model 1", ImVec2(-1, 0))) {

                viewer.data().clear();
                viewer.data().set_mesh(targetModel1.m_V, targetModel1.m_F);
                viewer.core.align_camera_center(targetModel1.m_V, targetModel1.m_F);
            }

            if (ImGui::Button("Show Control Cage 1", ImVec2(-1, 0))) {


                viewer.data().add_edges(controlCage1.headPoints, controlCage1.tailPoints, Eigen::RowVector3d(0, 1, 1));
            }

            // Add buttons
            if (ImGui::Button("Show Deformed Cage 1", ImVec2(-1, 0))) {
                viewer.data().add_edges(deformedCage1.headPoints, deformedCage1.tailPoints, Eigen::RowVector3d(1, 0.5, 0.5));
            }
            if (ImGui::Button("Deformed Model 1", ImVec2(-1, 0))) {
                viewer.data().clear();
                viewer.data().add_edges(deformedCage1.headPoints, deformedCage1.tailPoints, Eigen::RowVector3d(1,0.5,0.5));
                viewer.data().set_mesh(updatedVertex1, targetModel1.m_F);
                viewer.core.align_camera_center(updatedVertex1, targetModel1.m_F);
            }

            if (ImGui::Button("Collapse Deformed Cage 1", ImVec2(-1, 0))){
                viewer.data().clear();
                viewer.data().set_mesh(updatedVertex1, targetModel1.m_F);
                viewer.core.align_camera_center(updatedVertex1, targetModel1.m_F);
            }


        }

        //Add new group
        if (ImGui::CollapsingHeader("Deformation Model 2", ImGuiTreeNodeFlags_DefaultOpen)) {
            // Add buttons
            if (ImGui::Button("Show / Reset Model 2", ImVec2(-2, 0))) {
                viewer.data().clear();
                viewer.data().set_mesh(targetModel2.m_V, targetModel2.m_F);
                viewer.core.align_camera_center(targetModel2.m_V, targetModel2.m_F);
            }

            if (ImGui::Button("Show Control Cage 2", ImVec2(-2, 0))) {


                viewer.data().add_edges(controlCage2.headPoints, controlCage2.tailPoints, Eigen::RowVector3d(0, 1, 1));
            }

            // Add buttons
            if (ImGui::Button("Show Deformed Cage 2", ImVec2(-2, 0))) {
                viewer.data().add_edges(deformedCage2.headPoints, deformedCage2.tailPoints, Eigen::RowVector3d(1, 0.5, 0.5));
            }
            if (ImGui::Button("Deform Model 2", ImVec2(-2, 0))) {
                viewer.data().clear();
                viewer.data().add_edges(deformedCage2.headPoints, deformedCage2.tailPoints, Eigen::RowVector3d(1,0.5,0.5));
                viewer.data().set_mesh(updatedVertex2, targetModel2.m_F);
                viewer.core.align_camera_center(updatedVertex2, targetModel2.m_F);
            }

            if (ImGui::Button("Collapse Deformed Cage 2", ImVec2(-2, 0))){
                viewer.data().clear();
                viewer.data().set_mesh(updatedVertex2, targetModel2.m_F);
                viewer.core.align_camera_center(updatedVertex2, targetModel2.m_F);
            }


        }
        // Add new group
        if (ImGui::CollapsingHeader("Deformation Model 3", ImGuiTreeNodeFlags_DefaultOpen)) {
            // Add buttons
            if (ImGui::Button("Show / Reset Model 3", ImVec2(-3, 0))) {

                viewer.data().clear();
                viewer.data().set_mesh(targetModel3.m_V, targetModel3.m_F);
                viewer.core.align_camera_center(targetModel3.m_V, targetModel3.m_F);
            }

            if (ImGui::Button("Show Control Cage 3", ImVec2(-3, 0))) {


                viewer.data().add_edges(controlCage3.headPoints, controlCage3.tailPoints, Eigen::RowVector3d(0, 1, 1));
            }

            // Add buttons
            if (ImGui::Button("Show Deformed Cage 3", ImVec2(-3, 0))) {
                viewer.data().add_edges(deformedCage3.headPoints, deformedCage3.tailPoints, Eigen::RowVector3d(1, 0.5, 0.5));
            }
            if (ImGui::Button("Deform Model 3", ImVec2(-3, 0))) {
                viewer.data().clear();
                viewer.data().add_edges(deformedCage3.headPoints, deformedCage3.tailPoints, Eigen::RowVector3d(1,0.5,0.5));
                viewer.data().set_mesh(updatedVertex3, targetModel3.m_F);
                viewer.core.align_camera_center(updatedVertex3, targetModel3.m_F);
            }

            if (ImGui::Button("Collapse Deformed Cage 3", ImVec2(-3, 0))){
                viewer.data().clear();
                viewer.data().set_mesh(updatedVertex3, targetModel3.m_F);
                viewer.core.align_camera_center(updatedVertex3, targetModel3.m_F);
            }


        }
    };

    // registered a event handler
    viewer.callback_key_down = &key_down;
    viewer.launch();
}
