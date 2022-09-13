/*************************************************************************\

Sergei Grudinin, 2014
Mikhail Karasikov, 2014-2017
All Rights Reserved.

\**************************************************************************/


#include <string>
#include <vector>
#include <exception>
#include <cstdlib>
#include <tclap/CmdLine.h>

#include "engine.hpp"
#include "cClassicalParserPDB.hpp"
#include "cProteinMapper.hpp"
#include <iostream>
#include <fstream>
#include <chrono>


const double kDefaultAtomCutoffDistance = 7.0;
const double kDefaultResidueCutoffDistance = 5.0;
const double kDefaultHydrogenBondCutoffLength = 2.5;

const double kSaRadius = 2.0;  // 1.9 (+VdW!!!)
const double kSaSmoothedMargin = 1.0;  // length of weight change from 1 to 0 in Angstrom

#define RELEASE_VERSION 1

int main(int argc, char *argv[]) {
	std::string out_file;
	//std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
	try {
		TCLAP::CmdLine cmd("Protein Structure Prediction", ' ', "0.1");


		std::vector<std::string> allowed_modes = {
			"energy",
			"predict",
			"quality",
			"greedy",
			"featurize",
			"complete",
			"map",
			"gradmap",
            "sh"
		};
		TCLAP::ValuesConstraint<std::string> modes_constraint(allowed_modes);

		typedef TCLAP::ValueArg<std::string> CLField;

		CLField mode_arg("", "mode", "Mode", true, "", &modes_constraint, cmd);
		cmd.parse(std::min(argc, 3), argv);

		std::string mode = mode_arg.getValue();
		mode_arg.reset();

		cClassicalParserPDB parser;

        if (mode == "sh") {
            TCLAP::ValueArg<size_t> _order("p", "order", "Expansion order of spharical harmonics, 5 by default", false, 5, "int", cmd);
            CLField input_file_arg("i", "in_PDB", "PDB input file path", true, "", "path", cmd);
            CLField output_feature_file_arg("o", "out_maps", "Feature output file path", true, "", "path", cmd);
            TCLAP::SwitchArg skip_parsing_errors_arg("", "skip_errors", "Skip errors in PDB", cmd, false);

            cmd.parse(argc, argv);
            out_file = output_feature_file_arg.getValue();
            cProtein protein;

            if (!parser.parse(input_file_arg.getValue(), &protein,
//                              [](const cAtom &atom) { return atom.isHydrogen() || !atom.isBackbone(); },
                              [](const cAtom &atom) { return atom.isHydrogen(); },
                              skip_parsing_errors_arg.getValue(),
                              skip_parsing_errors_arg.getValue())
                ) return 1;
            protein.complete(false);

            if (!protein.isBackboneValid()) {
                std::cerr << "Error: Protein has incomplete backbone\n";
                exit(1);
            }

            if (!writeSHfeatures(&protein, out_file, _order.getValue())) {
                std::cerr << "Error when writing SH features file <" << out_file << ">\n";
                return 1;
            }


        }
		else if (mode == "map") {

			TCLAP::ValueArg<size_t> _gridSize("m", "grid_size", "Number of voxels in each dimension of the map", true, 0, "int", cmd);
			TCLAP::ValueArg<double> _atomSpreading("s", "spreading", "Spreading of the atoms or the residues in coarse mode (in A)", false, 1, "double", cmd);
			TCLAP::ValueArg<double> _voxWidth("v", "vox_size", "Size of the voxel (in A)", false, 1, "double", cmd);
			TCLAP::ValueArg<int> _types("t", "atom_types", "number of atom types", false, 1, "int", cmd);
			//TCLAP::ValueArg<bool> _hydrogen("h", "hydrogens", "consider hydrogen atoms", false, false, "bool", cmd);
			TCLAP::ValueArg<bool> _orient("", "orient", "orient the local maps", false, true, "bool", cmd);
			TCLAP::ValueArg<int> _neighbRes("", "skip_neighb", "skip the neighbour residues", false, 0, "int", cmd);

			CLField input_file_arg("i", "in_PDB", "PDB input file path", true, "", "path", cmd);
			TCLAP::SwitchArg skip_parsing_errors_arg("", "skip_errors", "Skip errors in PDB", cmd, false);
			TCLAP::SwitchArg native_arg("", "native", "Native PDB structure", cmd, false);
			TCLAP::SwitchArg coarse_arg("", "coarse", "coarse residue representation", cmd, false);
			CLField output_feature_file_arg("o", "out_maps", "Maps output file path", true, "", "path", cmd);
			cmd.parse(argc, argv);
			out_file = output_feature_file_arg.getValue();
			cProtein protein;

			if (!parser.parse(input_file_arg.getValue(), &protein,
				[](const cAtom &atom) { return atom.isHydrogen(); },
				skip_parsing_errors_arg.getValue(),
				skip_parsing_errors_arg.getValue()))
				return 1;
			protein.complete(false);

			std::unique_ptr<cProteinMapper> mapper;

			if (coarse_arg.getValue()) {
				mapper.reset(
					new cCoarseResidueMapper(_gridSize.getValue(),
						_atomSpreading.getValue(), _voxWidth.getValue(),  _orient.getValue())
				);
			}
			else {
				mapper.reset(
					new cAtomicResidueMapper(_types.getValue(),_gridSize.getValue(),
						_atomSpreading.getValue(), _voxWidth.getValue(), false, _orient.getValue(), _neighbRes.getValue())
				);
			}
			

			if (!protein.isBackboneValid()) {
				std::cerr << "Error: Protein has incomplete backbone\n";
				exit(1);
			}
			std::string scoreFileName = "";
			if (!native_arg.getValue()) {
				std::fstream f((input_file_arg.getValue()+".sco").c_str());
				if (f.good()) {
					scoreFileName = input_file_arg.getValue() + ".sco";
				}
			}

			//int ss_worked = system(("stride " + input_file_arg.getValue() + " > outss_txt").c_str());
			std::ifstream infile("outss_txt");
			std::string line;
			std::map<int, char> ssByRes;
			std::map<int, float> areaByRes;

			while (std::getline(infile, line))
			{
				if (line.size() > 26 &&
					line.at(0) == 'A' &&
					line.at(1) == 'S' &&
					line.at(2) == 'G') {
					int seqNumb = atoi(line.substr(10, 5).c_str());
					ssByRes.insert(std::pair<int, char>(seqNumb, line.at(24)));
					int area = atof(line.substr(63, 6).c_str());
					areaByRes.insert(std::pair<int, float>(seqNumb, area));
				}
			}
			if (!writeAllMaps(&protein, mapper.get(),
				output_feature_file_arg.getValue(), native_arg.getValue(), scoreFileName,
				ssByRes, areaByRes)) {
				std::cerr << "Error when writing maps file <" << output_feature_file_arg.getValue() << ">\n";
				return 1;
			}



		}
#if RELEASE_VERSION

		else {

			std::cerr << "Error, map is the only available mode for this version " << std::endl;
		}

#else//RELEASE_VERSION
        else if (mode == "gradmap") {

			TCLAP::ValueArg<size_t> _gridSize("m", "grid_size", "Number of voxels in each dimension of the map", true, 0, "int", cmd);
			TCLAP::ValueArg<double> _atomSpreading("s", "atom_spreading", "Spreading of the atoms (in A)", false, 1, "double", cmd);
			TCLAP::ValueArg<double> _voxWidth("v", "vox_size", "Size of the voxel (in A)", false, 1, "double", cmd);
			TCLAP::ValueArg<int> _types("t", "atom_types", "number of atom types", false, 1, "int", cmd);

			CLField input_file_arg("i", "in_PDB", "PDB input file path", true, "", "path", cmd);
			TCLAP::SwitchArg skip_parsing_errors_arg("", "skip_errors", "Skip errors in PDB", cmd, false);
			TCLAP::SwitchArg native_arg("", "native", "Native PDB structure", cmd, false);
			CLField gradient_map_file_arg("g", "grad_maps", "input gradient maps", true, "", "path", cmd);
			CLField output_grad_file_arg("o", "out_grad", "output gradient list", true, "", "path", cmd);
			cmd.parse(argc, argv);
			cProtein protein;

			if (!parser.parse(input_file_arg.getValue(), &protein,
				[](const cAtom &atom) { return atom.isHydrogen(); },
				skip_parsing_errors_arg.getValue(),
				skip_parsing_errors_arg.getValue()))
				return 1;

			std::unique_ptr<cAtomicResidueMapper> mapper;

			mapper.reset(
				new cAtomicResidueMapper(_gridSize.getValue(),
					_atomSpreading.getValue(), _voxWidth.getValue())
			);



			if (!protein.isBackboneValid()) {
				std::cerr << "Error: Protein has incomplete backbone\n";
				exit(1);
			}
			

			if (!computeGrad(&protein, mapper.get(),
				gradient_map_file_arg.getValue(), output_grad_file_arg.getValue())) {
				std::cerr << "Error whencomputing gradient file <" << output_grad_file_arg.getValue() << ">\n";
				return 1;
			}




		}
		else if (mode == "featurize") { // Features generation mode
			std::vector<std::string> allowed_featurization_modes = {
				"allatom", "backboneatom", "residues", "hbonds", "solvation"
			};
			TCLAP::ValuesConstraint<std::string> featurization_modes_constraint(allowed_featurization_modes);
			CLField featurization_mode_arg("", "featurization", "Featurization type", true, "", &featurization_modes_constraint, cmd);

			cmd.parse(std::min(argc, 5), argv);
			std::string featurization_mode = featurization_mode_arg.getValue();
			cmd.reset();

			cProtein protein;
			std::unique_ptr<cProteinFeaturizer> featurizer;

			CLField input_file_arg("i", "in_PDB", "PDB input file path", true, "", "path", cmd);
			TCLAP::SwitchArg skip_parsing_errors_arg("", "skip_errors", "Skip errors in PDB", cmd, false);
			CLField output_feature_file_arg("o", "out_features", "Features output file path", true, "", "path", cmd);

			if (featurization_mode == "residues") {
				TCLAP::ValueArg<size_t> num_DOF_arg("d", "num_DOF", "Number of parameters describing relative position of a residue pair", true, 0, "int", cmd);
				TCLAP::ValueArg<size_t> num_dist_bins_arg("b", "num_dist_bins", "Number of bins in histogram for distance", true, 0, "int", cmd);
				TCLAP::ValueArg<size_t> num_angle_bins_arg("a", "num_angle_bins", "Number of bins in histogram for each angle", true, 0, "int", cmd);
				TCLAP::ValueArg<double> cutoff_arg("c", "cutoff", "Cutoff for distance between two ends of side-chains", false, kDefaultResidueCutoffDistance, "double", cmd);
				TCLAP::ValueArg<size_t> skipped_neighbourhood_arg("n", "skipped_neighbourhood", "Skip neighbouring residues in topological distance up to |n|", false, 0, "int", cmd);
				TCLAP::ValueArg<double> smoothing_sigma_arg("s", "sigma", "Sigma for truncated gaussian smoothing (standard deviation in Angstroms)", false, 0.0, "double", cmd);
				cmd.parse(argc, argv);

				if (!parser.parse(input_file_arg.getValue(), &protein,
					[](const cAtom &atom) { return atom.isHydrogen() || !atom.isBackbone(); },
					skip_parsing_errors_arg.getValue(),
					skip_parsing_errors_arg.getValue()))
					return 1;
				protein.complete(true);

				featurizer.reset(
					new cResiduePairwiseHist(num_DOF_arg.getValue(),
						num_dist_bins_arg.getValue(), num_angle_bins_arg.getValue(),
						cutoff_arg.getValue(), skipped_neighbourhood_arg.getValue(),
						smoothing_sigma_arg.getValue())
				);
			}
			else if (featurization_mode == "allatom") {
				TCLAP::SwitchArg residue_type_dependent_arg("", "residue_type_dependent", "Grouping depends on the atom's residue type", cmd, false);
				TCLAP::SwitchArg atom_place_dependent_arg("", "atom_place_dependent", "Grouping depends on the atom's topological position in the residue", cmd, false);
				TCLAP::ValueArg<size_t> num_bins_arg("b", "num_bins", "Number of bins in histogram", true, 0, "int", cmd);
				TCLAP::ValueArg<double> cutoff_arg("c", "cutoff", "Cutoff for distance between a pair of atoms", false, kDefaultAtomCutoffDistance, "double", cmd);
				TCLAP::ValueArg<size_t> skipped_neighbourhood_arg("n", "skipped_neighbourhood", "Skip neighbouring atoms in topological distance up to |n|", false, 0, "int", cmd);
				TCLAP::ValueArg<double> smoothing_sigma_arg("s", "sigma", "Sigma for truncated gaussian smoothing (standard deviation in Angstroms)", false, 0.0, "double", cmd);
				cmd.parse(argc, argv);

				if (!parser.parse(input_file_arg.getValue(), &protein,
					[](const cAtom &atom) { return atom.isHydrogen(); },
					skip_parsing_errors_arg.getValue(),
					skip_parsing_errors_arg.getValue()))
					return 1;
				protein.complete(false);

				featurizer.reset(
					new cAtomPairwiseHist(residue_type_dependent_arg.getValue(),
						atom_place_dependent_arg.getValue(),
						num_bins_arg.getValue(),
						cutoff_arg.getValue(), skipped_neighbourhood_arg.getValue(),
						smoothing_sigma_arg.getValue())
				);
			}
			else if (featurization_mode == "backboneatom") {
				TCLAP::SwitchArg residue_type_dependent_arg("", "residue_type_dependent", "Grouping depends on the atom's residue type", cmd, false);
				TCLAP::SwitchArg atom_place_dependent_arg("", "atom_place_dependent", "Grouping depends on the atom's topological position in the residue", cmd, false);
				TCLAP::ValueArg<size_t> num_bins_arg("b", "num_bins", "Number of bins in histogram", true, 0, "int", cmd);
				TCLAP::ValueArg<double> cutoff_arg("c", "cutoff", "Cutoff for distance between a pair of atoms", false, kDefaultAtomCutoffDistance, "double", cmd);
				TCLAP::ValueArg<size_t> skipped_neighbourhood_arg("n", "skipped_neighbourhood", "Skip neighbouring atoms in topological distance up to |n|", false, 0, "int", cmd);
				TCLAP::ValueArg<double> smoothing_sigma_arg("s", "sigma", "Sigma for truncated gaussian smoothing (standard deviation in Angstroms)", false, 0.0, "double", cmd);
				cmd.parse(argc, argv);

				if (!parser.parse(input_file_arg.getValue(), &protein,
					[](const cAtom &atom) { return atom.isHydrogen() || !atom.isBackbone(); },
					skip_parsing_errors_arg.getValue(),
					skip_parsing_errors_arg.getValue()))
					return 1;
				protein.complete(false);

				featurizer.reset(
					new cAtomPairwiseHist(residue_type_dependent_arg.getValue(),
						atom_place_dependent_arg.getValue(),
						num_bins_arg.getValue(),
						cutoff_arg.getValue(), skipped_neighbourhood_arg.getValue(),
						smoothing_sigma_arg.getValue())
				);
			}
			else if (featurization_mode == "solvation") {
				TCLAP::ValueArg<size_t> num_dist_bins_arg("b", "num_dist_bins", "Number of bins in histogram for distance", true, 0, "int", cmd);
				TCLAP::ValueArg<size_t> num_angle_bins_arg("a", "num_angle_bins", "Number of bins in histogram for angle", true, 0, "int", cmd);
				TCLAP::ValueArg<double> cutoff_arg("c", "cutoff", "Cutoff for distance between residue and solvent", false, kDefaultAtomCutoffDistance, "double", cmd);
				TCLAP::ValueArg<double> smoothing_sigma_arg("s", "sigma", "Sigma for truncated gaussian smoothing (standard deviation in Angstroms)", false, 0.0, "double", cmd);
				cmd.parse(argc, argv);

				if (!parser.parse(input_file_arg.getValue(), &protein,
					[](const cAtom &atom) { return atom.isHydrogen() || !atom.isBackbone(); },
					skip_parsing_errors_arg.getValue(),
					skip_parsing_errors_arg.getValue()))
					return 1;
				protein.complete(true);

				featurizer.reset(
					new cSolvationShellHist(kSaRadius,
						kSaSmoothedMargin,
						num_dist_bins_arg.getValue(),
						num_angle_bins_arg.getValue(),
						cutoff_arg.getValue(),
						smoothing_sigma_arg.getValue())
				);
			}
			else if (featurization_mode == "hbonds") {
				TCLAP::ValueArg<size_t> num_dist_bins_arg("b", "num_dist_bins", "Number of bins in histogram for distance", true, 0, "int", cmd);
				TCLAP::ValueArg<size_t> num_angle_bins_arg("a", "num_angle_bins", "Number of bins in histogram for each angle", true, 0, "int", cmd);
				TCLAP::ValueArg<double> cutoff_arg("c", "cutoff", "Cutoff for length of hydrogen bond", false, kDefaultHydrogenBondCutoffLength, "double", cmd);
				TCLAP::ValueArg<size_t> skipped_neighbourhood_arg("n", "skipped_neighbourhood", "Skip hydrogen bonds for neighbouring residues in topological distance up to |n|", false, 0, "int", cmd);
				TCLAP::ValueArg<double> smoothing_sigma_arg("s", "sigma", "Sigma for truncated gaussian smoothing (standard deviation in Angstroms)", false, 0.0, "double", cmd);
				cmd.parse(argc, argv);

				if (!parser.parse(input_file_arg.getValue(), &protein,
					[](const cAtom &atom) { return atom.isHydrogen() || !atom.isBackbone(); },
					skip_parsing_errors_arg.getValue(),
					skip_parsing_errors_arg.getValue()))
					return 1;
				protein.complete(true);

				featurizer.reset(
					new cHydrogenBondsHist(num_dist_bins_arg.getValue(), num_angle_bins_arg.getValue(),
						cutoff_arg.getValue(), skipped_neighbourhood_arg.getValue(),
						smoothing_sigma_arg.getValue())
				);
			}
			else {
				assert(false);
			}

			if (!protein.isBackboneValid()) {
				std::cerr << "Error: Protein has incomplete backbone\n";
				exit(1);
			}

			if (!writeFeatureVector(&protein, featurizer.get(),
				output_feature_file_arg.getValue())) {
				std::cerr << "Error when writing features file <" << output_feature_file_arg.getValue() << ">\n";
				return 1;
			}
		}
		else if (mode == "energy") { // Energy matrix generation mode
			CLField input_file("i", "in_PDB", "PDB input file path", true, "", "path", cmd);
			CLField rotamer_lib_name("r", "rotLib", "Rotamer Library path", true, "", "path", cmd);
			CLField output_mat_file("m", "outMAT", "MAT output file path", true, "", "path", cmd);

			cmd.parse(argc, argv);

			cProtein protein;
			if (!parser.parse(input_file.getValue(), &protein))
				return 1;
			if (!protein.initRotamersFromDunbrackBDLib(rotamer_lib_name.getValue()))
				return 1;

			protein.complete(true); // with Hydrogens

			if (!writeEnergyMatrix(&protein, output_mat_file.getValue())) {
				std::cerr << "Error when writing mat file <" << output_mat_file.getValue() << ">\n";
				return 1;
			}
		}
		else if (mode == "predict") {

			CLField input_file("i", "in_PDB", "PDB input file path", true, "", "path", cmd);
			CLField rotamer_lib_name("r", "rotLib", "Rotamer Library path", true, "", "path", cmd);
			CLField rotamers_file("", "rotamers", "Path to ASCII input file with rotamers", true, "", "path", cmd);
			CLField output_file("o", "outPDB", "PDB output file path", true, "", "path", cmd);

			cmd.parse(argc, argv);

			cProtein protein;
			if (!parser.parse(input_file.getValue(), &protein))
				return 1;
			if (!protein.initRotamersFromDunbrackBDLib(rotamer_lib_name.getValue()))
				return 1;

			protein.complete(); // without hydrogens

			if (!protein.setRotamers(rotamers_file.getValue()))
				return 1;

			if (!protein.writePDB(output_file.getValue()))
				return 1;

		}
		else if (mode == "quality") { // test quality

			CLField testing_file("i", "in_PDB", "Path to testing PDB input file", true, "", "path", cmd);
			CLField ethalon_file("e", "ethalonPDB", "Path to ethalon PDB input file", true, "", "path", cmd);

			cmd.parse(argc, argv);

			cProtein testing;
			if (!parser.parse(testing_file.getValue(), &testing))
				exit(1);
			testing.complete();

			cProtein ethalon;
			if (!parser.parse(ethalon_file.getValue(), &ethalon))
				exit(1);
			ethalon.complete();

			printf("Chi_1 quality: %3.1lf%%\n", chi1Quality(testing, ethalon));
			printf("Chi_1+2 quality: %3.1lf%%\n", chi12Quality(testing, ethalon));
			printf("RMSD: %lf\n", rmsdQuality(testing, ethalon));
			printf("Clashes: %zd from %zd\n", clashQuality(testing), clashQuality(ethalon));
			printf("Energy: %lf from %lf\n", energyQuality(&testing), energyQuality(&ethalon));

		}
		else if (mode == "greedy") {

			CLField input_file("i", "in_PDB", "PDB input file path", true, "", "path", cmd);
			CLField rotamer_lib_name("r", "rotLib", "Rotamer Library path", true, "", "path", cmd);
			CLField output_file("o", "outPDB", "PDB output file path", true, "", "path", cmd);

			cmd.parse(argc, argv);

			cProtein protein;
			if (!parser.parse(input_file.getValue(), &protein))
				return 1;
			if (!protein.initRotamersFromDunbrackBDLib(rotamer_lib_name.getValue()))
				return 1;

			protein.complete();

			greedySearch(&protein, 10);

			protein.writePDB(output_file.getValue());

		}
		else if (mode == "complete") {

			CLField input_file_arg("i", "in_PDB", "PDB input file path", true, "", "path", cmd);
			CLField output_file_arg("o", "outPDB", "PDB output file path", true, "", "path", cmd);

			cmd.parse(argc, argv);

			cProtein protein;
			if (!parser.parse(input_file_arg.getValue(), &protein, [](const cAtom &atom ) { return atom.isHydrogen(); } ))
				return 1;

			protein.complete(true);
			protein.writePDB(output_file_arg.getValue(), true);
		}
#endif //RELEASE_VERSION
	}
	catch (TCLAP::ArgException &e) {
		std::cerr << "ERROR: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}
	catch (std::exception &e) {
		std::cerr << "ERROR: " << e.what() << std::endl;
		std::remove(out_file.c_str());
		return 1;
	}
	catch (...) {
		std::cerr << "Unknown exception" << std::endl;
		return 1;
	}

	//std::cout << "Total time : " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() -start).count() << " ms" <<std::endl;
	return 0;
}
