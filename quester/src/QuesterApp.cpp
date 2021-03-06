/**
 *
 *
 */

#include "QuesterApp.hpp"


namespace quester {

QuesterApp::QuesterApp (int argc, char** argv) {
    parse_args (argc, argv);
}

void QuesterApp::run() {
    if (opts.print_help_message == true) {
        print_usage ();
    } else {

    }
}


void QuesterApp::parse_args (int argc, char** argv) {
    namespace po = boost::program_options;
    namespace fs = boost::filesystem;

    std::string config_file;
    std::string verbosity;
    std::string quiteness;

    bool batch_mode = false;
    bool interactive_mode = false;
    std::string script_name;
    std::string command;

    exe_name = fs::path(std::string(argv[0])).stem().string();

    po::options_description options("Options");
    options.add_options()
        ("help,h",        
                "Produce help message")
        ("config,c", po::value <std::string> (&config_file),
                "Config file name")
        ("version",
                "Print version string")
        ("run-command,r", po::value<std::string>(&command), 
                "Command to run ... e.g., run filename.py")
        ("script_file,f", po::value<std::string>(&script_name), 
                "Script file to run.")
        ("batch,b", po::bool_switch(&batch_mode), 
                "Run in batch mode.")
        ("interactive,i", po::bool_switch(&interactive_mode), 
                "Run in interactive mode -- much like MATLAB command mode.")
        ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
        ("long-version", 
                "Print the version info in long, verbose format.")
        ;        

    // command line options
    po::options_description cmd_line_opts;
    cmd_line_opts.add (options).add (hidden_options);

    // config file options
    po::options_description conf_file_opts;
    conf_file_opts.add (options).add (hidden_options);

    // shown in the help message
    po::options_description usage (80, 40);
    usage.add (options);
    std::stringstream out;
    out << "Usage: " << exe_name << " [options] " << std::endl;
    out << usage;
    usage_message = out.str();

    // do the actual parsing
    try {
        po::variables_map users_choices;
        po::store (po::command_line_parser (argc, argv)
                        .options (cmd_line_opts).run(), users_choices);
        po::notify (users_choices);

        if (config_file.length() > 0) {
            std::ifstream file (config_file);
            po::store (po::parse_config_file(file, usage), users_choices);
        }
        po::notify (users_choices);

        // perform option checking
        if (batch_mode && interactive_mode) {
            usage_error_message = "Interactive and batch mode are mutually exclusive, please chose one of them.";
            opts.print_help_message = true;
            return;
        }

        if (script_name.length () > 0 && command.length () > 0) {
            usage_error_message = "Options --run-command and --script_file are mutually exclusive, please choose one of them.";
            opts.print_help_message = true;
            return;
        }


        if (script_name.length () > 0) {
            opts.command = "run \"" + script_name + "\"";
        }

        if (command.length () > 0) {
            opts.command = command;
        }

        if (batch_mode && opts.command.length () == 0) {
            usage_error_message = "Batch mode was requested however, neither command nor script name was specified, aborting.";
            opts.print_help_message = true;
            return;
        }

        if (!batch_mode && !interactive_mode) {
            if (script_name.length () > 0) {
                opts.mode = Options::Mode::Batch;
            } else if (command.length () > 0) {
                opts.mode = Options::Mode::Interactive;
            }
        } else if (batch_mode) {
            opts.mode = Options::Mode::Batch;
        } else if (interactive_mode) {
            opts.mode = Options::Mode::Interactive;
        }
        
        if (users_choices.count ("help")) {
            opts.print_help_message = true;
        }


    } catch (const po::unknown_option & e) {
        usage_error_message = e.what();
        opts.print_help_message = true;
    }
}

void QuesterApp::print_usage () {
    if (usage_error_message.length () > 0) {
        std::cout << std::endl << "ERROR: " << usage_error_message 
                  << std::endl << std::endl;
    }

    //std::cout << get_banner() << endl;
    std::cout << usage_message << std::endl;
}

}

