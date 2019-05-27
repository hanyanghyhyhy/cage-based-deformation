// DF: Deformation/ Deformation Cage
// CF: Control Point/ Control Vetex
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/adjacency_list.h>
#include <igl/edge_lengths.h>
#include <igl/jet.h>
#include <Eigen/SparseLU>
#include "nanoflann.hpp"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/edges.h>
#include <igl/file_exists.h>

double EPSILON = 1.0e-19;
Eigen::MatrixXd headPointsCage, trailPointsCage, headPointsCageDF, trailPointsCageDF;
Eigen::MatrixXd mvcWeights1;

Eigen::MatrixXd operator*(Eigen::MatrixXd weights, Eigen::MatrixXd cageVetex);

inline double computeDeterminant(Eigen::MatrixXd u0, Eigen::MatrixXd u1, Eigen::MatrixXd u2) {
    double determinant = u0(0) * (u1(1) * u2(2) - u2(1) * u1(2))
                         - u0(1) * (u1(0) * u2(2) - u2(0) * u1(2))
                         + u0(2) * (u1(0) * u2(1) - u2(0) * u1(1));
    return determinant;
}

class TargetModel {
public:
    TargetModel(std::string fileName) : nv_len(0.5), point_size(8), line_width(0.5), sel_vidx(7), mode(0) {
        //initial vertices and faces
        if (!igl::file_exists(fileName)) {
            std::cout << "[error] cannot locate model file at " << fileName << "\nPress any key to exit\n";
            char c;
            std::cin >> c;
            exit(1);
        }

        igl::readOBJ(fileName, m_V, m_F);
        // calculate VN
        igl::per_face_normals(m_V, m_F, m_FN);
    }

    ~TargetModel() {}

    Eigen::MatrixXd m_V; // mode 1 -p
    Eigen::MatrixXi m_F;
    Eigen::MatrixXd m_FN;
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

    void showModel(igl::opengl::glfw::Viewer &viewer) {

        viewer.data().set_mesh(m_V, m_F);
    }

private:
};

class ControlCage {
public:
    ControlCage(std::string fileName) {
        if (!igl::file_exists(fileName)) {
            std::cout << "[error] cannot locate model file at " << fileName << "\nPress any key to exit\n";
            char c;
            std::cin >> c;
            exit(1);
        }

        igl::readOBJ(fileName, m_V, m_F);
    }

    ~ControlCage() {}

    Eigen::MatrixXd m_V;
    Eigen::MatrixXi m_F;

    void loadCage(Eigen::MatrixXd &headPoints, Eigen::MatrixXd &trailPoints) {
        // The edges function returns the head and the trail position of the edges
        // use the position vectors to retrieve the edges of the cage
        // store them in headPoints and trailPoints.

        Eigen::MatrixXi Edges;
        igl::edges(m_F, Edges);
        headPoints = Eigen::MatrixXd(Edges.rows(), 3);
        trailPoints = Eigen::MatrixXd(Edges.rows(), 3);

        for (int i = 0; i < Edges.rows(); i++) {
            headPoints.row(i) = m_V.row(Edges(i, 0));
            trailPoints.row(i) = m_V.row(Edges(i, 1));
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
Eigen::MatrixXd
computeBarycentricCoordinates(int index, int rows, Eigen::MatrixXi controlFace, double theta0, double theta1,
                              double theta2, double l0, double l1, double l2) {

    Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(rows, 1);
    weights(controlFace(index, 0)) = sin(theta0) * l1 * l2;
    weights(controlFace(index, 1)) = sin(theta1) * l0 * l2;
    weights(controlFace(index, 2)) = sin(theta2) * l0 * l1;

    return weights;
}

// This function only compute the weight between single model point and control point.
Eigen::MatrixXd
computeSingleWeight(Eigen::MatrixXd modelVetex, Eigen::MatrixXd controlVetex, Eigen::MatrixXi controlFace) {

    Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(1, controlVetex.rows());

    // checks the distance between control point and model point
    for (int i = 0; i < controlVetex.rows(); i++) {

        Eigen::MatrixXd individualCP = controlVetex.row(i);

        Eigen::Vector3d d(modelVetex(0) - individualCP(0), modelVetex(1) - individualCP(1),
                          modelVetex(2) - individualCP(2));

        // if the distance between control vetex and model vetex is too close, return 1 to weight.
        if (d.norm() <= EPSILON) {
            weights(i) = 1.0;
            return weights;
        }
    }

    for (int j = 0; j < controlFace.rows(); j++) {
//        double w0,w1,w2;

        Eigen::MatrixXd p0 = controlVetex.row(controlFace(j, 0));
        Eigen::MatrixXd p1 = controlVetex.row(controlFace(j, 1));
        Eigen::MatrixXd p2 = controlVetex.row(controlFace(j, 2));

        Eigen::Vector3d d0(modelVetex(0) - p0(0), modelVetex(1) - p0(1), modelVetex(2) - p0(2));
        Eigen::Vector3d d1(modelVetex(0) - p1(0), modelVetex(1) - p1(1), modelVetex(2) - p1(2));
        Eigen::Vector3d d2(modelVetex(0) - p2(0), modelVetex(1) - p2(1), modelVetex(2) - p2(2));

        Eigen::Vector3d u0 = d0.normalized();
        Eigen::Vector3d u1 = d1.normalized();
        Eigen::Vector3d u2 = d2.normalized();

        double l0 = (u1 - u2).norm();
        double l1 = (u0 - u2).norm();
        double l2 = (u0 - u1).norm();

        double theta0 = 2 * asin(l0 * 0.5);
        double theta1 = 2 * asin(l1 * 0.5);
        double theta2 = 2 * asin(l2 * 0.5);

        double h = (theta0 + theta1 + theta2) * 0.5;

        if ((M_PI - h) <= EPSILON) {
            // x lie on t, use 2D barycentric coordinates
            return computeBarycentricCoordinates(j, controlVetex.rows(), controlFace, theta0, theta1, theta2, l0, l1,
                                                 l2);
        }


        double c0 = (2 * sin(h) * sin(h - theta0)) / (sin(theta1) * sin(theta2)) - 1.0;
        double c1 = (2 * sin(h) * sin(h - theta1)) / (sin(theta0) * sin(theta2)) - 1.0;
        double c2 = (2 * sin(h) * sin(h - theta2)) / (sin(theta0) * sin(theta1)) - 1.0;

        double determinant = computeDeterminant(u0, u1, u2);

        double sign = (determinant < EPSILON) ? -1 : 1;

        double s0 = sign * determinant * sqrt(1.0 - c0 * c0);
        double s1 = sign * determinant * sqrt(1.0 - c1 * c1);
        double s2 = sign * determinant * sqrt(1.0 - c2 * c2);

        // x lies outside t on the same plane, ignore t
        if (s0 <= EPSILON || s1 <= EPSILON || s2 <= EPSILON) {
            continue;
        }

        weights(controlFace(j, 0)) +=
                (theta0 - c1 * theta2 - c2 * theta1) / (2 * d0.norm() * sin(theta2) * sqrt(1 - c1 * c1));
        weights(controlFace(j, 1)) +=
                (theta1 - c0 * theta2 - c2 * theta0) / (2 * d1.norm() * sin(theta0) * sqrt(1 - c2 * c2));
        weights(controlFace(j, 2)) +=
                (theta2 - c1 * theta0 - c0 * theta1) / (2 * d2.norm() * sin(theta1) * sqrt(1 - c0 * c0));

    }

    return weights;
}

// Compute the weights between model veteice and cage before deformation
// This function could be called before the deformation to expedite the computation
// For each vetex on the model surface, we compute the weighs between it and all the control points
Eigen::MatrixXd computeWeight(Eigen::MatrixXd modelVetex, Eigen::MatrixXd cageVetex, Eigen::MatrixXi cageFace) {
    Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(modelVetex.rows(), cageVetex.rows());
    for (int i = 0; i < modelVetex.rows(); i++) {
        weights.row(i) = computeSingleWeight(modelVetex.row(i), cageVetex, cageFace);
    }
    return weights;
}

// After deformation, update the vetices in model using computed weights.
Eigen::MatrixXd updateDeformedVetex(Eigen::MatrixXd weights, Eigen::MatrixXd cageVetex) {
    Eigen::MatrixXd updatedVetex = Eigen::MatrixXd::Zero(weights.rows(), 3);

    for (int i = 0; i < weights.rows(); i++) {
        Eigen::MatrixXd singleRow = Eigen::MatrixXd::Zero(1, 3);
        double weightSum = 0.0;
        for (int j = 0; j < weights.cols(); j++) {
            singleRow += weights(i, j) * cageVetex.row(j);
            weightSum += weights(i, j);
        }
        updatedVetex.row(i) = singleRow / weightSum;
    }

    return updatedVetex;
}

// Use operator overload
Eigen::MatrixXd operator*(Eigen::MatrixXd weights, Eigen::MatrixXd cageVetex) {
    Eigen::MatrixXd updatedVetex = Eigen::MatrixXd::Zero(weights.rows(), 3);

    for (int i = 0; i < weights.rows(); i++) {
        Eigen::MatrixXd singleRow = Eigen::MatrixXd::Zero(1, 3);
        double weightSum = 0.0;
        for (int j = 0; j < weights.cols(); j++) {
            singleRow += weights(i, j) * cageVetex.row(j);
            weightSum += weights(i, j);
        }
        updatedVetex.row(i) = singleRow / weightSum;
    }

    return updatedVetex;
}

int main(int argc, char *argv[]) {

    std::string modelPath1 = "../meshes/Beast.obj";
    std::string cagePath1 = "../meshes/Beast_Cage.obj";
    std::string dfPath1 = "../meshes/Beast_Cage_Deformed.obj";
    std::string modePath2 = "../meshes/Bench.obj";
    std::string cagePath2 = "../meshes/Bench_Cage.obj";
    std::string dfPath2 = "../meshes/Bench_Cage_deformed.obj";
    std::string modelPath3 = "../meshes/cactus.obj";
    std::string cagePath3 = "../meshes/cactus_cage.onj";

    // load the model
    TargetModel targetModel1 = TargetModel(modelPath1);
    ControlCage controlCage1 = ControlCage(cagePath1);
    // Calculate the weight in advance to exoedite the computation
    mvcWeights1 = computeWeight(targetModel1.m_V, controlCage1.m_V, controlCage1.m_F);





    // Init the viewer
    igl::opengl::glfw::Viewer viewer;
    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    menu.callback_draw_viewer_menu = [&]() {
        // Add new group
        if (ImGui::CollapsingHeader("Deformation 1", ImGuiTreeNodeFlags_DefaultOpen)) {
            // Add buttons
            if (ImGui::Button("Show/Reset Model", ImVec2(-1, 0))) {

                viewer.data().clear();

                viewer.data().set_mesh(targetModel1.m_V, targetModel1.m_F);
                viewer.core.align_camera_center(targetModel1.m_V, targetModel1.m_F);
            }

            if (ImGui::Button("Show Control Cage", ImVec2(-1, 0))) {

                controlCage1.loadCage(headPointsCage, trailPointsCage);
                viewer.data().add_edges(headPointsCage, trailPointsCage, Eigen::RowVector3d(0, 1, 1));
            }

            // Add buttons
            if (ImGui::Button("Show Deformed Cage", ImVec2(-1, 0))) {
                ControlCage deformedCage1 = ControlCage(dfPath1);
                deformedCage1.loadCage(headPointsCageDF, trailPointsCageDF);
                viewer.data().add_edges(headPointsCageDF, trailPointsCageDF, Eigen::RowVector3d(1, 0.5, 0.5));
            }
            if (ImGui::Button("Deform Model", ImVec2(-1, 0))) {
                ControlCage deformedCage1 = ControlCage(dfPath1);
                deformedCage1.loadCage(headPointsCageDF, trailPointsCageDF);
                Eigen::MatrixXd updatedVertex = mvcWeights1 * deformedCage1.m_V;
                viewer.data().set_mesh(updatedVertex, targetModel1.m_F);
                viewer.core.align_camera_center(updatedVertex, targetModel1.m_F);
            }


        }

    };

    // registered a event handler
    viewer.callback_key_down = &key_down;
    viewer.launch();


}
