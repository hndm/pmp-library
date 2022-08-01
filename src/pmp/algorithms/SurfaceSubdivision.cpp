// Copyright 2011-2022 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "pmp/algorithms/SurfaceSubdivision.h"
#include "pmp/algorithms/DifferentialGeometry.h"

namespace pmp {

SurfaceSubdivision::SurfaceSubdivision(SurfaceMesh& mesh) : mesh_(mesh)
{
    points_ = mesh_.vertex_property<Point>("v:point");
    vfeature_ = mesh_.get_vertex_property<bool>("v:feature");
    efeature_ = mesh_.get_edge_property<bool>("e:feature");
}

void SurfaceSubdivision::catmull_clark()
{
    // reserve memory
    size_t nv = mesh_.n_vertices();
    size_t ne = mesh_.n_edges();
    size_t nf = mesh_.n_faces();
    mesh_.reserve(nv + ne + nf, 2 * ne + 4 * nf, 4 * nf);

    // get properties
    auto vpoint = mesh_.add_vertex_property<Point>("catmull:vpoint");
    auto epoint = mesh_.add_edge_property<Point>("catmull:epoint");
    auto fpoint = mesh_.add_face_property<Point>("catmull:fpoint");

    // compute face vertices
    for (auto f : mesh_.faces())
    {
        fpoint[f] = centroid(mesh_, f);
    }

    // compute edge vertices
    for (auto e : mesh_.edges())
    {
        // boundary or feature edge?
        if (mesh_.is_boundary(e) || (efeature_ && efeature_[e]))
        {
            epoint[e] = 0.5f * (points_[mesh_.vertex(e, 0)] +
                                points_[mesh_.vertex(e, 1)]);
        }

        // interior edge
        else
        {
            Point p(0, 0, 0);
            p += points_[mesh_.vertex(e, 0)];
            p += points_[mesh_.vertex(e, 1)];
            p += fpoint[mesh_.face(e, 0)];
            p += fpoint[mesh_.face(e, 1)];
            p *= 0.25f;
            epoint[e] = p;
        }
    }

    // compute new positions for old vertices
    for (auto v : mesh_.vertices())
    {
        // isolated vertex?
        if (mesh_.is_isolated(v))
        {
            vpoint[v] = points_[v];
        }

        // boundary vertex?
        else if (mesh_.is_boundary(v))
        {
            auto h1 = mesh_.halfedge(v);
            auto h0 = mesh_.prev_halfedge(h1);

            Point p = points_[v];
            p *= 6.0;
            p += points_[mesh_.to_vertex(h1)];
            p += points_[mesh_.from_vertex(h0)];
            p *= 0.125;

            vpoint[v] = p;
        }

        // interior feature vertex?
        else if (vfeature_ && vfeature_[v])
        {
            Point p = points_[v];
            p *= 6.0;
            int count(0);

            for (auto h : mesh_.halfedges(v))
            {
                if (efeature_[mesh_.edge(h)])
                {
                    p += points_[mesh_.to_vertex(h)];
                    ++count;
                }
            }

            if (count == 2) // vertex is on feature edge
            {
                p *= 0.125;
                vpoint[v] = p;
            }
            else // keep fixed
            {
                vpoint[v] = points_[v];
            }
        }

        // interior vertex
        else
        {
            // weights from SIGGRAPH paper "Subdivision Surfaces in Character Animation"

            const Scalar k = mesh_.valence(v);
            Point p(0, 0, 0);

            for (auto vv : mesh_.vertices(v))
                p += points_[vv];

            for (auto f : mesh_.faces(v))
                p += fpoint[f];

            p /= (k * k);

            p += ((k - 2.0f) / k) * points_[v];

            vpoint[v] = p;
        }
    }

    // assign new positions to old vertices
    for (auto v : mesh_.vertices())
    {
        points_[v] = vpoint[v];
    }

    // split edges
    for (auto e : mesh_.edges())
    {
        // feature edge?
        if (efeature_ && efeature_[e])
        {
            auto h = mesh_.insert_vertex(e, epoint[e]);
            auto v = mesh_.to_vertex(h);
            auto e0 = mesh_.edge(h);
            auto e1 = mesh_.edge(mesh_.next_halfedge(h));

            vfeature_[v] = true;
            efeature_[e0] = true;
            efeature_[e1] = true;
        }

        // normal edge
        else
        {
            mesh_.insert_vertex(e, epoint[e]);
        }
    }

    // split faces
    for (auto f : mesh_.faces())
    {
        auto h0 = mesh_.halfedge(f);
        mesh_.insert_edge(h0, mesh_.next_halfedge(mesh_.next_halfedge(h0)));

        auto h1 = mesh_.next_halfedge(h0);
        mesh_.insert_vertex(mesh_.edge(h1), fpoint[f]);

        auto h =
            mesh_.next_halfedge(mesh_.next_halfedge(mesh_.next_halfedge(h1)));
        while (h != h0)
        {
            mesh_.insert_edge(h1, h);
            h = mesh_.next_halfedge(
                mesh_.next_halfedge(mesh_.next_halfedge(h1)));
        }
    }

    // clean-up properties
    mesh_.remove_vertex_property(vpoint);
    mesh_.remove_edge_property(epoint);
    mesh_.remove_face_property(fpoint);
}

void SurfaceSubdivision::loop()
{
    if (!mesh_.is_triangle_mesh())
    {
        auto what = "SurfaceSubdivision: Not a triangle mesh.";
        throw InvalidInputException(what);
    }

    // reserve memory
    size_t nv = mesh_.n_vertices();
    size_t ne = mesh_.n_edges();
    size_t nf = mesh_.n_faces();
    mesh_.reserve(nv + ne, 2 * ne + 3 * nf, 4 * nf);

    // add properties
    auto vpoint = mesh_.add_vertex_property<Point>("loop:vpoint");
    auto epoint = mesh_.add_edge_property<Point>("loop:epoint");

    // compute vertex positions
    for (auto v : mesh_.vertices())
    {
        // isolated vertex?
        if (mesh_.is_isolated(v))
        {
            vpoint[v] = points_[v];
        }

        // boundary vertex?
        else if (mesh_.is_boundary(v))
        {
            auto h1 = mesh_.halfedge(v);
            auto h0 = mesh_.prev_halfedge(h1);

            Point p = points_[v];
            p *= 6.0;
            p += points_[mesh_.to_vertex(h1)];
            p += points_[mesh_.from_vertex(h0)];
            p *= 0.125;
            vpoint[v] = p;
        }

        // interior feature vertex?
        else if (vfeature_ && vfeature_[v])
        {
            Point p = points_[v];
            p *= 6.0;
            int count(0);

            for (auto h : mesh_.halfedges(v))
            {
                if (efeature_[mesh_.edge(h)])
                {
                    p += points_[mesh_.to_vertex(h)];
                    ++count;
                }
            }

            if (count == 2) // vertex is on feature edge
            {
                p *= 0.125;
                vpoint[v] = p;
            }
            else // keep fixed
            {
                vpoint[v] = points_[v];
            }
        }

        // interior vertex
        else
        {
            Point p(0, 0, 0);
            Scalar k(0);

            for (auto vv : mesh_.vertices(v))
            {
                p += points_[vv];
                ++k;
            }
            p /= k;

            Scalar beta =
                (0.625 - pow(0.375 + 0.25 * std::cos(2.0 * M_PI / k), 2.0));

            vpoint[v] = points_[v] * (Scalar)(1.0 - beta) + beta * p;
        }
    }

    // compute edge positions
    for (auto e : mesh_.edges())
    {
        // boundary or feature edge?
        if (mesh_.is_boundary(e) || (efeature_ && efeature_[e]))
        {
            epoint[e] =
                (points_[mesh_.vertex(e, 0)] + points_[mesh_.vertex(e, 1)]) *
                Scalar(0.5);
        }

        // interior edge
        else
        {
            auto h0 = mesh_.halfedge(e, 0);
            auto h1 = mesh_.halfedge(e, 1);
            Point p = points_[mesh_.to_vertex(h0)];
            p += points_[mesh_.to_vertex(h1)];
            p *= 3.0;
            p += points_[mesh_.to_vertex(mesh_.next_halfedge(h0))];
            p += points_[mesh_.to_vertex(mesh_.next_halfedge(h1))];
            p *= 0.125;
            epoint[e] = p;
        }
    }

    // set new vertex positions
    for (auto v : mesh_.vertices())
    {
        points_[v] = vpoint[v];
    }

    // insert new vertices on edges
    for (auto e : mesh_.edges())
    {
        // feature edge?
        if (efeature_ && efeature_[e])
        {
            auto h = mesh_.insert_vertex(e, epoint[e]);
            auto v = mesh_.to_vertex(h);
            auto e0 = mesh_.edge(h);
            auto e1 = mesh_.edge(mesh_.next_halfedge(h));

            vfeature_[v] = true;
            efeature_[e0] = true;
            efeature_[e1] = true;
        }

        // normal edge
        else
        {
            mesh_.insert_vertex(e, epoint[e]);
        }
    }

    // split faces
    Halfedge h;
    for (auto f : mesh_.faces())
    {
        h = mesh_.halfedge(f);
        mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
        h = mesh_.next_halfedge(h);
        mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
        h = mesh_.next_halfedge(h);
        mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
    }

    // clean-up properties
    mesh_.remove_vertex_property(vpoint);
    mesh_.remove_edge_property(epoint);
}

void SurfaceSubdivision::quad_tri()
{
    // split each edge evenly into two parts
    for (auto e : mesh_.edges())
    {
        mesh_.insert_vertex(e, 0.5f * (points_[mesh_.vertex(e, 0)] +
                                       points_[mesh_.vertex(e, 1)]));
    }

    // subdivide faces without repositioning
    for (auto f : mesh_.faces())
    {
        size_t f_val = mesh_.valence(f) / 2;
        if (f_val == 3)
        {
            // face was a triangle
            Halfedge h0 = mesh_.halfedge(f);
            Halfedge h1 = mesh_.next_halfedge(mesh_.next_halfedge(h0));
            mesh_.insert_edge(h0, h1);

            h0 = mesh_.next_halfedge(h0);
            h1 = mesh_.next_halfedge(mesh_.next_halfedge(h0));
            mesh_.insert_edge(h0, h1);

            h0 = mesh_.next_halfedge(h0);
            h1 = mesh_.next_halfedge(mesh_.next_halfedge(h0));
            mesh_.insert_edge(h0, h1);
        }
        else
        {
            // quadrangulate the rest
            Halfedge h0 = mesh_.halfedge(f);
            Halfedge h1 = mesh_.next_halfedge(mesh_.next_halfedge(h0));
            h1 = mesh_.insert_edge(h0, h1);
            mesh_.insert_vertex(mesh_.edge(h1), centroid(mesh_, f));

            auto h = mesh_.next_halfedge(
                mesh_.next_halfedge(mesh_.next_halfedge(h1)));
            while (h != h0)
            {
                mesh_.insert_edge(h1, h);
                h = mesh_.next_halfedge(
                    mesh_.next_halfedge(mesh_.next_halfedge(h1)));
            }
        }
    }

    auto new_pos =
        mesh_.add_vertex_property<Point>("quad_tri:new_position", Point(0));

    for (auto v : mesh_.vertices())
    {
        if (mesh_.is_boundary(v))
        {
            new_pos[v] = 0.5 * points_[v];

            // add neighbouring vertices on boundary
            for (auto vv : mesh_.vertices(v))
            {
                if (mesh_.is_boundary(vv))
                {
                    new_pos[v] += 0.25 * points_[vv];
                }
            }
        }
        else
        {
            // count the number of faces and quads surrounding the vertex
            int n_faces = 0, n_quads = 0;
            for (auto f : mesh_.faces(v))
            {
                n_faces++;
                n_quads += mesh_.valence(f) == 4;
            }

            if (n_quads == 0)
            {
                // vertex is surrounded only by triangles
                double a =
                    2.0 * pow(3.0 / 8.0 +
                                  (std::cos(2.0 * M_PI / n_faces) - 1.0) / 4.0,
                              2.0);
                double b = (1.0 - a) / n_faces;

                new_pos[v] = a * points_[v];
                for (auto vv : mesh_.vertices(v))
                {
                    new_pos[v] += b * points_[vv];
                }
            }
            else if (n_quads == n_faces)
            {
                // vertex is surrounded only by quads
                double c = (n_faces - 3.0) / n_faces;
                double d = 2.0 / pow(n_faces, 2.0);
                double e = 1.0 / pow(n_faces, 2.0);

                new_pos[v] = c * points_[v];
                for (auto h : mesh_.halfedges(v))
                {
                    new_pos[v] += d * points_[mesh_.to_vertex(h)];
                    new_pos[v] +=
                        e * points_[mesh_.to_vertex(mesh_.next_halfedge(h))];
                }
            }
            else
            {
                // vertex is surrounded by triangles and quads
                double alpha = 1.0 / (1.0 + 0.5 * n_faces + 0.25 * n_quads);
                double beta = 0.5 * alpha;
                double gamma = 0.25 * alpha;

                new_pos[v] = alpha * points_[v];
                for (auto h : mesh_.halfedges(v))
                {
                    new_pos[v] += beta * points_[mesh_.to_vertex(h)];
                    if (mesh_.valence(mesh_.face(h)) == 4)
                    {
                        new_pos[v] +=
                            gamma *
                            points_[mesh_.to_vertex(mesh_.next_halfedge(h))];
                    }
                }
            }
        }
    }

    // apply new positions to the mesh
    for (auto v : mesh_.vertices())
    {
        points_[v] = new_pos[v];
    }

    mesh_.remove_vertex_property(new_pos);
}

void SurfaceSubdivision::quad_tri_interpolating(int n_steps, int n_iterations)
{
    auto vsharp_edges_ = mesh_.vertex_property<int>("v:sharp_edges", 0);
    auto esharpness_ = mesh_.edge_property<Scalar>("e:sharpness", -1.0);

    size_t nv = mesh_.n_vertices();

    // define interpolating points
    Eigen::MatrixX3d ipoints, cpoints;
    ipoints.resize(nv, Eigen::NoChange);
    int row_idx = 0;
    for (Vertex v : mesh_.vertices())
    {
        ipoints.row(row_idx) = static_cast<Eigen::Vector3d>(points_[v]);
        ++row_idx;
    }
    cpoints = ipoints;

    // subdivision matrix
    Eigen::SparseMatrix<double> W(nv, nv);
    W.setIdentity();

    Eigen::SparseMatrix<double> P, V;
    std::vector<Eigen::Triplet<double>> P_trip, V_trip;
    for (int i = 0; i < n_steps; ++i)
    {
        P_trip.clear();
        V_trip.clear();

        const size_t n_initial_vertices = mesh_.n_vertices();

        auto vidx = mesh_.add_vertex_property<IndexType>("quad_tri:vidx");
        int idx = 0;

        // vertex point weights for linear subdivision matrix
        for (Vertex v : mesh_.vertices())
        {
            P_trip.emplace_back(idx, idx, 1.0);
            vidx[v] = idx++;
        }

        // edge point weights for linear subdivsion matrix
        for (Edge e : mesh_.edges())
        {
            P_trip.emplace_back(idx, vidx[mesh_.vertex(e, 0)], 0.5);
            P_trip.emplace_back(idx, vidx[mesh_.vertex(e, 1)], 0.5);
            ++idx;
        }

        // face point weights for linear subdivision matrix
        for (Face f : mesh_.faces())
        {
            size_t fval = mesh_.valence(f);
            if (fval >= 4)
            {
                for (Vertex v : mesh_.vertices(f))
                {
                    P_trip.emplace_back(idx, vidx[v], 1.0 / fval);
                }
                ++idx;
            }
        }

        idx = n_initial_vertices;

        // split edges linearly
        for (auto e : mesh_.edges())
        {
            // for linear subdivision
            auto epoint =
                0.5 * (points_[mesh_.vertex(e, 0)] + points_[mesh_.vertex(e, 1)]);
            Halfedge h = mesh_.insert_vertex(e, epoint);

            vidx[mesh_.to_vertex(h)] = idx++;

            vsharp_edges_[mesh_.to_vertex(h)] = 2;
            esharpness_[mesh_.edge(h)] = esharpness_[e];
            esharpness_[mesh_.edge(mesh_.next_halfedge(h))] = esharpness_[e];
        }


        // split faces linearly
        for (Face f : mesh_.faces())
        {
            size_t fVal = mesh_.valence(f) / 2; // initial valence of face
            if (fVal == 3)
            {
                auto h = mesh_.halfedge(f);
                mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
                h = mesh_.next_halfedge(h);
                mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
                h = mesh_.next_halfedge(h);
                mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
            }
            else if (fVal >= 4)
            {
                auto fpoint = centroid(mesh_, f);

                auto h0 = mesh_.halfedge(f);
                mesh_.insert_edge(h0, mesh_.next_halfedge(mesh_.next_halfedge(h0)));

                auto h1 = mesh_.next_halfedge(h0);
                mesh_.insert_vertex(mesh_.edge(h1), fpoint);

                vidx[mesh_.to_vertex(h1)] = idx++;

                auto h = mesh_.next_halfedge(
                    mesh_.next_halfedge(mesh_.next_halfedge(h1)));
                while (h != h0)
                {
                    mesh_.insert_edge(h1, h);
                    h = mesh_.next_halfedge(
                        mesh_.next_halfedge(mesh_.next_halfedge(h1)));
                }
            }
        }

        // compute smoothing mask for linearly subdivided mesh
        for (Vertex v : mesh_.vertices())
        {
            // isolated vertex
            if (mesh_.is_isolated(v))
            {
                V_trip.emplace_back(vidx[v], vidx[v], 1.0);
            }

            // boundary vertex
            else if (mesh_.is_boundary(v))
            {
                V_trip.emplace_back(vidx[v], vidx[v], 0.5);
                for (Vertex vv : mesh_.vertices(v))
                {
                    if (mesh_.is_boundary(vv))
                    {
                        V_trip.emplace_back(vidx[v], vidx[vv], 0.25);
                    }
                }
            }

            // corner vertex
            else if (vsharp_edges_[v] > 2)
            {
                V_trip.emplace_back(vidx[v], vidx[v], 1.0);
            }
            // smooth / sharp / semi-sharp vertex
            else
            {
                std::vector<Eigen::Triplet<double>> smooth_weights, sharp_weights;
                float vsharpness = 0.0f;

                // compute weights for sharp vertex
                if (vsharp_edges_[v] == 2)
                {
                    sharp_weights.emplace_back(vidx[v], vidx[v], 0.5);

                    for (auto h : mesh_.halfedges(v))
                    {
                        if (esharpness_[mesh_.edge(h)] > 0.0f)
                        {
                            sharp_weights.emplace_back(
                                vidx[v], vidx[mesh_.to_vertex(h)], 0.25);
                            vsharpness += esharpness_[mesh_.edge(h)];
                        }
                    }

                    vsharpness *= 0.5f;
                }

                // crease vertex
                if (vsharpness >= 1.0f)
                {
                    V_trip.insert(V_trip.end(), sharp_weights.begin(),
                                  sharp_weights.end());
                }
                // compute weights for smooth vertex
                else
                {
                    int n_edges = 0, n_quads = 0;
                    for (Face f : mesh_.faces(v))
                    {
                        n_edges++;
                        n_quads += mesh_.valence(f) == 4;
                    }

                    // surrounded by triangles
                    if (n_quads == 0)
                    {
                        double a =
                            2.0 * pow(3.0 / 8.0 +
                                          (std::cos(2.0 * M_PI / n_edges)) / 4.0,
                                      2.0) -
                            0.25;
                        double b = (1.0 - a) / n_edges;

                        smooth_weights.emplace_back(vidx[v], vidx[v], a);
                        for (Vertex vv : mesh_.vertices(v))
                        {
                            smooth_weights.emplace_back(vidx[v], vidx[vv], b);
                        }
                    }
                    // surrounded by quads
                    else if (n_quads == n_edges)
                    {
                        double c = (n_edges - 3.0) / n_edges;
                        double d = 2.0 / pow(n_edges, 2.0);
                        double e = 1.0 / pow(n_edges, 2.0);

                        smooth_weights.emplace_back(vidx[v], vidx[v], c);
                        for (Halfedge h : mesh_.halfedges(v))
                        {
                            smooth_weights.emplace_back(
                                vidx[v], vidx[mesh_.to_vertex(h)], d);
                            smooth_weights.emplace_back(
                                vidx[v],
                                vidx[mesh_.to_vertex(mesh_.next_halfedge(h))], e);
                        }
                    }
                    // surrounded by triangles and quads
                    else
                    {
                        double alpha = 1.0 / (1.0 + n_edges / 2.0 + n_quads / 4.0);
                        double beta = alpha / 2.0;
                        double gamma = alpha / 4.0;

                        smooth_weights.emplace_back(vidx[v], vidx[v], alpha);
                        for (Halfedge h : mesh_.halfedges(v))
                        {
                            smooth_weights.emplace_back(
                                vidx[v], vidx[mesh_.to_vertex(h)], beta);
                            if (mesh_.valence(mesh_.face(h)) == 4)
                            {
                                smooth_weights.emplace_back(
                                    vidx[v],
                                    vidx[mesh_.to_vertex(mesh_.next_halfedge(h))],
                                    gamma);
                            }
                        }
                    }

                    // semi-sharp vertex
                    if (vsharpness > 0.0f)
                    {
                        for (auto t : sharp_weights)
                            V_trip.emplace_back(t.row(), t.col(),
                                                t.value() * vsharpness);
                        for (auto t : smooth_weights)
                            V_trip.emplace_back(t.row(), t.col(),
                                                t.value() * (1.0 - vsharpness));
                    }
                    // smooth vertex
                    else
                    {
                        V_trip.insert(V_trip.end(), smooth_weights.begin(),
                                      smooth_weights.end());
                    }
                }
            }
        }

        // update edge sharpness
        for (Edge e : mesh_.edges())
        {
            if (esharpness_[e] > 0.0)
            {
                esharpness_[e] = std::max(esharpness_[e] - 1.0, 0.0);
                if (esharpness_[e] <= 0.0)
                {
                    vsharp_edges_[mesh_.vertex(e, 0)] -= 1;
                    vsharp_edges_[mesh_.vertex(e, 1)] -= 1;
                }
            }
        }

        mesh_.remove_vertex_property(vidx);

        P.resize(mesh_.n_vertices(), n_initial_vertices);
        P.setFromTriplets(P_trip.begin(), P_trip.end());

        V.resize(mesh_.n_vertices(), mesh_.n_vertices());
        V.setFromTriplets(V_trip.begin(), V_trip.end());

        W = V * P * W;
    }

    // update positions for initial vertices
    for (int k = 0; k < n_iterations; ++k)
    {
        cpoints += ipoints - (W * cpoints).topRows(nv);
    }

    std::cout << ipoints - (W * cpoints).topRows(nv) << std::endl;

    // compute vertex positions for subdivided mesh
    Eigen::MatrixX3d spoints = W * cpoints;

    row_idx = 0;
    for (Vertex v : mesh_.vertices())
    {
        points_[v] = spoints.row(row_idx);
        ++row_idx;
    }
}
} // namespace pmp
