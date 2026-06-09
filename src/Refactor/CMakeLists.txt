cmake_minimum_required(VERSION 4.0)
project(Refactor)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

add_executable(Refactor
        main.cpp
        Clusters/ClusterMethod.h
        Clusters/Opticluster.cpp
        Clusters/Opticluster.h
        DataExporters/ClusterExport.cpp
        DataExporters/ClusterExport.h
        DataStructures/ListVector.cpp
        DataStructures/ListVector.cpp
        DataStructures/ListVector.h
        DataStructures/RAbundVector.cpp
        DataStructures/RAbundVector.h
        MothurDependencies/Utils.cpp
        MothurDependencies/Utils.h
        DataExporters/ClusterData.cpp
        DataExporters/ClusterData.h
        DataExporters/OptiClusterData.cpp
        DataExporters/OptiClusterData.h
        Clusters/Cluster.cpp
        Clusters/Cluster.h
        DataStructures/SparseDistanceMatrix.cpp
        DataStructures/SparseDistanceMatrix.h
        DataStructures/PDistCell.h
        DataStructures/PDistCellMin.h
        DataStructures/OptiData.cpp
        DataStructures/OptiData.h
        Clusters/ClusterMetric.h
        Clusters/Metrics/accuracy.h
        Clusters/Metrics/f1score.h
        Clusters/Metrics/fdr.h
        Clusters/Metrics/fn.h
        Clusters/Metrics/fp.h
        Clusters/Metrics/fpfn.h
        Clusters/Metrics/mcc.h
        Clusters/Metrics/npv.h
        Clusters/Metrics/ppv.h
        Clusters/Metrics/sensitivity.h
        Clusters/Metrics/specificity.h
        Clusters/Metrics/tn.h
        Clusters/Metrics/tp.h
        Clusters/Metrics/tptn.h)
