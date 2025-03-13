cmake_minimum_required(VERSION 3.10)
project(CeresTest)

# 設置 C++ 標準
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 尋找需要的包
find_package(Ceres REQUIRED)
# glog 可能沒有 CMake 配置文件，使用 PkgConfig 或直接鏈接
find_package(glog QUIET)
if(NOT glog_FOUND)
    # 嘗試使用 pkg-config
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(GLOG libglog QUIET)
    endif()
endif()
find_package(absl REQUIRED)

# 添加可執行文件
add_executable(ceres_test test/ceres_test.cpp)

# 連結庫
target_link_libraries(ceres_test 
    Ceres::ceres 
    # 根據 glog 查找結果使用不同的鏈接方式
    $<$<BOOL:${glog_FOUND}>:glog::glog>
    $<$<BOOL:${GLOG_FOUND}>:${GLOG_LIBRARIES}>
    $<$<NOT:$<OR:$<BOOL:${glog_FOUND}>,$<BOOL:${GLOG_FOUND}>>>:glog>
    absl::base
    absl::strings
    absl::log_internal_message
)

# 設置包含目錄
target_include_directories(ceres_test PRIVATE 
    ${CERES_INCLUDE_DIRS}
)

# 添加可執行文件
add_executable(10bar src/ceres_for_10bar.cpp)

# 連結庫
target_link_libraries(10bar 
    Ceres::ceres 
    # 根據 glog 查找結果使用不同的鏈接方式
    $<$<BOOL:${glog_FOUND}>:glog::glog>
    $<$<BOOL:${GLOG_FOUND}>:${GLOG_LIBRARIES}>
    $<$<NOT:$<OR:$<BOOL:${glog_FOUND}>,$<BOOL:${GLOG_FOUND}>>>:glog>
    absl::base
    absl::strings
    absl::log_internal_message
)

# 設置包含目錄
target_include_directories(10bar PRIVATE 
    ${CERES_INCLUDE_DIRS}
    include
)


# 輸出編譯信息
message(STATUS "CERES_INCLUDE_DIRS: ${CERES_INCLUDE_DIRS}")
message(STATUS "CERES_LIBRARIES: ${CERES_LIBRARIES}")