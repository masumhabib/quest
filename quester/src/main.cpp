/* 
 * File:   main.cpp
 * @brief: QUEST main program that can run simulations in interactive
 *         or batch mode.
 *
 * Copyright (C) 2016 K M Masum Habib <masum.habib@ail.com>
 *
 */

#include "QuesterApp.hpp"
//#include "read_python.h"
//#include <quest.hpp>
#include <iostream>
//#include <vector>



//using namespace utils::stds;
//using utils::strings::ttos;
//using utils::python::create_main_module;
//using utils::python::read_option;
//using utils::python::write_option;
//using utils::python::eval_options;

//namespace bp = boost::python;


/*
 *  Entry point for the QUEST main program. 
 */
int main(int argc, char** argv) {
    using namespace quester;
    QuesterApp app (argc, argv);
    try {
        app.run();
    } catch (...) {
        //FIXME: print error message
        std::cout << "Unknown error occured" << std::endl;
    }

/*    wall_clock timer;
    timer.tic();

    string simu_file_name;
    if (argc < 2){
        cout << endl << " ERROR: Missing simulation file name.";
        cout << endl << " Usage: tmf simulation_file.py" << endl;
        exit(-1);
    }else{
        simu_file_name = argv[1];
    }    
    
    Workers workers;
    
    setVerbosity(10);
    vout.myId(workers.MyId());
    dout.myId(workers.MyId());
        
    string msg = greet();
    vout << msg;
            
    try{
        
        // Read python file for simulation parameters.
        vout << endl << " Reading simulation file " << simu_file_name << " ... ";
        Py_Initialize();
        bp::object main_module;
        create_main_module(main_module);
        write_option<bool>(main_module, "AmIMaster", workers.AmIMaster());
        eval_options(simu_file_name, main_module);
        vout << "done.";

        // split path, extract file names
        string out_file_name = simu_file_name.substr(0, simu_file_name.rfind(".py")) + ".dat";
        read_option<string>(main_module, "out_file_name", out_file_name);
        string::size_type pathpos = out_file_name.rfind("/");
        string out_path;
        string out_file_ext;
        if(pathpos == string::npos){
            out_path = "./";
        }else{
            out_path = out_file_name.substr(0, pathpos);
            out_file_name = out_file_name.substr(pathpos);
        }
        pathpos = out_file_name.rfind(".");
        if(pathpos == string::npos){
            out_file_ext = ".dat";
        }else{
            out_file_ext = out_file_name.substr(pathpos);
            out_file_name = out_file_name.substr(0, pathpos);
        }
        
        vector<double> Bz = {0.0};
        read_option<double>(main_module, "Bz", Bz);

        // set lattice parameter and k.p fermion parameter.
        GrapheneKpParams kpp;
        double a;
        if(read_option<double>(main_module, "a", a)){
            kpp.a(a);
        }
        bool bvalue;
        if(read_option<bool>(main_module, "is_neg_K", bvalue)){
            if(bvalue){
                kpp.K(-kpp.K());
            }
        } 

        // create hall_bar object
        hall_bar hb(workers, kpp);

        // output settings.
        if (read_option<bool>(main_module, "Taa_enabled", bvalue)){
            hb.enable_Taa(bvalue);
        }
        if (read_option<bool>(main_module, "Tab_enabled", bvalue)){
            hb.enable_Tab(bvalue);
        }
        if (read_option<bool>(main_module, "Tac_enabled", bvalue)){
            hb.enable_Tac(bvalue);
        }
        if (read_option<bool>(main_module, "Tad_enabled", bvalue)){
            hb.enable_Tad(bvalue);
        }
        if (read_option<bool>(main_module, "Tcb_enabled", bvalue)){
            hb.enable_Tcb(bvalue);
        }
        if (read_option<bool>(main_module, "Tcd_enabled", bvalue)){
            hb.enable_Tcd(bvalue);
        }
        if (read_option<bool>(main_module, "Tbd_enabled", bvalue)){
            hb.enable_Tbd(bvalue);
        }

        // energy grid
        double Emin = -1, Emax = 1, dE = 0.05;
        read_option<double>(main_module, "Emin", Emin);
        read_option<double>(main_module, "Emax", Emax);
        read_option<double>(main_module, "dE", dE);
        hb.set_energy_grid(Emin, Emax, dE);

        // geometry
        vout << endl << " Creating geometry      ... ";
        double length, width, pos;
        if(read_option<double>(main_module, "length", length) 
                && read_option<double>(main_module, "width", width)){
            hb.add_device(length, width);
        }else{
            throw runtime_error(" In main(): need to specify length and width of device.");
        }

        if (read_option<double>(main_module, "width_contact_a", width) 
                && read_option<double>(main_module, "position_contact_a", pos)){
            hb.add_contact_a(width, pos);
        }
        if (read_option<double>(main_module, "width_contact_b", width) 
                && read_option<double>(main_module, "position_contact_b", pos)){
            hb.add_contact_b(width, pos);
        }
        if (read_option<double>(main_module, "width_contact_c", width) 
                && read_option<double>(main_module, "position_contact_c", pos)){
            hb.add_contact_c(width, pos);
        }
        if (read_option<double>(main_module, "width_contact_d", width) 
                && read_option<double>(main_module, "position_contact_d", pos)){
            hb.add_contact_d(width, pos);
        }
        
        bool export_geometry = false;
        read_option<bool>(main_module, "export_geometry", export_geometry);
        if (export_geometry){
            hb.export_geometry(out_path + "/"
                + simu_file_name.substr(0, simu_file_name.rfind(".py")) + ".gjf");
        }

        vout << "done.";

        // gates
        double split;
        if(read_option<double>(main_module, "split", split)){
            vout << endl << " Adding gates ... ";
            hb.add_gate(split);
            vout << "done.";
        }

        // gate roughness
        // simple roughness from rms value
        double roughness_rms = 0;
        if(read_option<double>(main_module, "roughness_rms", roughness_rms)){
            vout << endl << " Adding roughness ... ";
            vout << endl << "   RMS = " << roughness_rms << " angstroms" << endl;
            hb.add_gate_roughness(roughness_rms);
            vout << "done.";
        }

        // roughness given as vectors
        vector<double> roughness_left_edge;
        vector<double> roughness_right_edge;
        read_option<double>(main_module, "roughness_left_edge", roughness_left_edge);
        read_option<double>(main_module, "roughness_right_edge", roughness_right_edge);
        if(!roughness_left_edge.empty() || !roughness_right_edge.empty()){
            vout << endl << " Adding roughness ... ";
            hb.add_gate_roughness(roughness_left_edge, roughness_right_edge);
            vout << "done.";
        }

        vout << endl << " Setting potential      ... ";   
        double VGl = 0, VGr = 0;
        read_option<double>(main_module, "VGl", VGl);
        read_option<double>(main_module, "VGr", VGr);
        hb.set_potential(VGl, VGr);
        bool export_potential = false;
        read_option<bool>(main_module, "export_potential", export_potential);
        if (export_potential){
            hb.export_potential(out_path + "/"
                + simu_file_name.substr(0, simu_file_name.rfind(".py")) + ".pot");
        }
        vout << "done.";

        vout << endl << " Simulation parameters  ... ";    
        vout << hb; 
        vout << endl << " ------------------------------------------------------------------";

        for (auto B:Bz){
            vout << endl << endl;
            kpp.Bz(B);

            vout << endl << " Setting Bz = " << B << " T";    
            vout << endl << " Generating Hamiltonian ... ";    
            hb.generate_ham();
            vout << "done.";

            vout << endl << " Running simulation on " << workers.N() << " CPU(s) ... " << endl;  
            bool dry_run = false;
            read_option<bool>(main_module, "dry_run", dry_run);

            if(!dry_run){
                hb.run();
                //hb.test();
            }
            vout << endl << "  done.";

            vout << endl << " Saving data to ";   
            stringstream out_file_name_full;
            out_file_name_full << out_file_name << "_B" << std::setprecision(3) << std::fixed << B;
            out_file_name_full << out_file_ext;

            vout << out_path << out_file_name_full.str() << " ... ";
            if(!dry_run){
                hb.save(out_path + "/" + out_file_name_full.str());
            }
            vout << "done.";
        }
        
        // save a copy of python script.
        std::ifstream    in_file(simu_file_name);
        std::ofstream    out_file(out_path + "/" + simu_file_name.substr(0, simu_file_name.rfind(".py")) + ".bkp.py");
        out_file << in_file.rdbuf();
        
    }catch(bp::error_already_set &e){
        dout << endl << "ERROR: In file " << simu_file_name << ": " << endl;
        PyErr_Print();
    }catch(std::exception &e){
        dout << endl << "ERROR: " << e.what() << endl;
    }
    
    vout << endl << " Runtime: " << ttos(timer.toc()) << endl;
    
    */
    return 0;
}

