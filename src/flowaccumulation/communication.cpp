#include "communication.h"
#include <paradem/tool.h>

#include <assert.h>
#include <fstream>
#include <iostream>

#define _unused(x) ((void)x)

void CommISend(msg_type& msg, int dest, int tag) {
    MPI_Request request;
    bytes_sent += msg.size();
    MPI_Isend(msg.data(), msg.size(), MPI_BYTE, dest, tag, MPI_COMM_WORLD, &request);
}

void CommInit(int* argc, char*** argv) {
    MPI_Init(argc, argv);
}

int CommRank() {
    int rank;
    int ret = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    assert(ret == MPI_SUCCESS);
    _unused(ret);
    return rank;
}

int CommSize() {
    int size;
    int ret = MPI_Comm_size(MPI_COMM_WORLD, &size);
    assert(ret == MPI_SUCCESS);
    _unused(ret);
    return size;
}

comm_count_type CommBytesSent() {
    return bytes_sent;
}
comm_count_type CommBytesRecv() {
    return bytes_recv;
}
void CommBytesReset() {
    bytes_recv = 0;
    bytes_sent = 0;
}

int CommGetTag(int from) {
    MPI_Status status;
    MPI_Probe(from, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    return status.MPI_TAG;
}

void Preparer(const GridInfo& gridDEMInfo, const std::vector<TileInfo>& tileDEMInfos, const GridInfo& gridDirInfo, const std::vector<TileInfo>& tileDirInfos, ObjectFactory* pIObjFactory,
              std::string filename) {

    Timer timer_overall;
    timer_overall.start();

    int jobs_out = 0;

    std::vector<msg_type> msgs;
    int size = CommSize();
    const int active_consumer_limit = size - 1;  //当前通信数
    int good_to_go = 1;
    CommBroadcast(&good_to_go, 0);

    for (size_t i = 0; i < tileDirInfos.size(); i++) {
        msgs.push_back(CommPrepare(&gridDEMInfo, &gridDirInfo, &tileDEMInfos[i], &tileDirInfos[i], &filename));
        CommISend(msgs.back(), (i % active_consumer_limit) + 1, TagFirst);
        jobs_out++;
    }
    std::cout << "all jobs out is " << jobs_out << std::endl;
    Grid<std::shared_ptr<IConsumer2Producer>> gridIConsumer2Producer;
    gridIConsumer2Producer.init(gridDirInfo.gridHeight, gridDirInfo.gridWidth);
    while (jobs_out--) {
        Consumer2Producer pC2P;
        CommRecv(&pC2P, nullptr, -1);
        gridIConsumer2Producer.at(pC2P.gridRow, pC2P.gridCol) = std::make_shared<Consumer2Producer>(pC2P);
    }

    std::cerr << "n First stage Tx = " << CommBytesSent() << " B" << std::endl;
    std::cerr << "n First stage Rx = " << CommBytesRecv() << " B" << std::endl;
    CommBytesReset();

    TimeInfo time_first_total;
    for (int row = 0; row < gridDEMInfo.gridHeight; row++)
        for (int col = 0; col < gridDEMInfo.gridWidth; col++)
            time_first_total += ((Consumer2Producer*)gridIConsumer2Producer.at(row, col).get())->time_info;
    std::shared_ptr<IProducer> pIProducer = pIObjFactory->createProducer();
    pIProducer->process(gridDirInfo, tileDirInfos, gridIConsumer2Producer);
    jobs_out = 0;
    for (size_t i = 0; i < tileDirInfos.size(); i++) {
        TileInfo tileDirInfo = tileDirInfos[i];
        std::shared_ptr<IProducer2Consumer> p2c = pIProducer->toConsumer(gridIConsumer2Producer.at(tileDirInfo.gridRow, tileDirInfo.gridCol).get());
        Producer2Consumer* pP2c = (Producer2Consumer*)p2c.get();
        CommSend(&gridDEMInfo, &gridDirInfo, &tileDEMInfos[i], &tileDirInfos[i], &filename, pP2c, (i % active_consumer_limit) + 1, TagSecond);
        jobs_out++;
    }
    TimeInfo time_second_total;
    while (jobs_out--) {
        TimeInfo temp;
        CommRecv(&temp, nullptr, -1);
        time_second_total += temp;
    }

    for (int i = 1; i < size; i++) {
        int temp;
        CommSend(&temp, nullptr, i, TagNull);
    }

    timer_overall.stop();

    std::cerr << "t First stage total overall time = " << time_first_total.overall << " s" << std::endl;
    std::cerr << "t First stage total io time = " << time_first_total.io << " s" << std::endl;
    std::cerr << "t First stage total calc time = " << time_first_total.calc << " s" << std::endl;
    std::cerr << "r First stage peak child VmPeak = " << time_first_total.vmpeak << std::endl;
    std::cerr << "r First stage peak child VmHWM = " << time_first_total.vmhwm << std::endl;

    std::cerr << "n Second stage Tx = " << CommBytesSent() << " B" << std::endl;
    std::cerr << "n Second stage Rx = " << CommBytesRecv() << " B" << std::endl;

    std::cerr << "t Second stage total overall time = " << time_second_total.overall << " s" << std::endl;
    std::cerr << "t Second stage total IO time = " << time_second_total.io << " s" << std::endl;
    std::cerr << "t Second stage total calc time = " << time_second_total.calc << " s" << std::endl;
    std::cerr << "r Second stage peak child VmPeak = " << time_second_total.vmpeak << std::endl;
    std::cerr << "r Second stage peak child VmHWM = " << time_second_total.vmhwm << std::endl;

    std::cerr << "t Producer overall time = " << timer_overall.elapsed() << " s" << std::endl;
    std::cerr << "t Producer calc time = " << ((Producer*)pIProducer.get())->timer_calc.elapsed() << " s" << std::endl;
    long vmpeak, vmhwm;
    ProcessMemUsage(vmpeak, vmhwm);
    std::cerr << "r Producer's VmPeak = " << vmpeak << std::endl;
    std::cerr << "r Producer's VmHWM = " << vmhwm << std::endl;
    std::string txtPath = gridDirInfo.outputFolder + "\\" + "gridInfo.txt";
    std::ofstream fout;
    fout.open(txtPath, std::ofstream::app);
    if (fout.fail()) {
        std::cout << "Open " << txtPath << " error!" << std::endl;
        return;
    }
    fout << gridDirInfo.tileHeight << std::endl;
    fout << gridDirInfo.tileWidth << std::endl;
    fout << gridDirInfo.gridHeight << std::endl;
    fout << gridDirInfo.gridWidth << std::endl;
    fout << gridDirInfo.grandHeight << std::endl;
    fout << gridDirInfo.grandWidth << std::endl;
    fout << gridDirInfo.cellSize << std::endl;
    fout.close();

    ////merge tiles
    // std::string path = gridDEMInfo.outputFolder;
    // std::string pathAccum = path + "\\gridInfo.txt";
    // GridInfo gridAccumInfo;
    // readGridInfo(pathAccum.data(), gridAccumInfo);
    // gridAccumInfo.inputFolder = path;
    // gridAccumInfo.outputFolder = path;
    // mergeTiles(gridAccumInfo);
}