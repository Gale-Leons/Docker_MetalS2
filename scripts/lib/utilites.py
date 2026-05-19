#!/usr/bin/env python

import getopt
import math
import os
import sys

sys.path.append(os.path.join(os.path.expanduser("~"), "development"))
from lib.extractSites import getSitesFromPdbFile


def processCommandLineArguments(argv):
    input1, input2, output, maxDist, rmPdb, scoreOnly = parseCommandLineOptions(argv)
    sites1, sites2 = handleInputs(input1, input2)
    pathToOutput = handleOutputs(output)

    return sites1, sites2, pathToOutput, maxDist, rmPdb, scoreOnly


def handleInputs(input1, input2):
    type1 = input1.get("type")
    type2 = input2.get("type")

    path1 = input1.get("path")
    path2 = input2.get("path")

    metalID1 = input1.get("metal")
    metalID2 = input2.get("metal")

    if not opened(path1) or not opened(path2):
        print("\nFiles do not exist or not accessible\n")
        sys.exit()

    sites1 = []
    sites2 = []

    # Both files are pdb
    if type1 == "pdb" and type2 == "pdb":
        sites1 = getSitesFromPdbFile(path1, metalID1)
        sites2 = getSitesFromPdbFile(path2, metalID2)

    else:
        print("\nFiles are in a bad format.\n")
        sys.exit()

    return sites1, sites2


def handleOutputs(output):
    if not output:
        outputDir = "results"
        outputPath = os.path.join(os.getcwd(), outputDir)
    else:
        try:
            output.split("/")[-2]
        except IndexError:
            outputPath = os.path.join(os.getcwd(), output)
        else:
            outputPath = output

    return outputPath


def parseCommandLineOptions(argv):
    # parse command line options
    try:
        opts, args = getopt.getopt(
            argv[1:],
            "d:hu",
            [
                "help",
                "usage",
                "qp=",
                "tp=",
                "qs=",
                "ts=",
                "qm=",
                "tm=",
                "rm_pdb",
                "score_only",
            ],
        )

    except getopt.GetoptError as err:
        print(
            "\n" + str(err)
        )  # will print something like "option -a is not recognized"

        if err.opt == "pq":
            print("Did you mean --qp?")
        if err.opt == "pt":
            print("Did you mean --tp?")
        if err.opt == "sq":
            print("Did you mean --qs?")
        if err.opt == "st":
            print("Did you mean --ts?")

        if err.opt == "mq":
            print("Did you mean --qm?")
        if err.opt == "mt":
            print("Did you mean --tm?")

        print("\nFor help use -h or --help\n")
        sys.exit()

    if opts:
        inputType1 = None
        inputType2 = None
        metalID1 = None
        metalID2 = None
        maxDist = 2.0
        rmPdb = False
        scoreOnly = False

        # process options
        for o, a in opts:
            if o == "--qp":
                inputType1 = "pdb"
                pathToInput1 = a
            elif o == "--tp":
                inputType2 = "pdb"
                pathToInput2 = a
            elif o == "--qs":
                inputType1 = "site"
                pathToInput1 = a
            elif o == "--ts":
                inputType2 = "site"
                pathToInput2 = a

            elif o == "--qm":
                metalID1 = int(a)
            elif o == "--tm":
                metalID2 = int(a)

            elif o == "--rm_pdb":
                rmPdb = True
            elif o == "--score_only":
                scoreOnly = True

            elif o == "-d":
                maxDist = math.fabs(float(a))
                if maxDist == 0:
                    print("Distance value must be higher than zero.")
                    sys.exit()

            elif o in ("-u", "--usage"):
                # print usage information and exit
                usage()
                sys.exit()
            elif o in ("-h", "--help"):
                # print help information and exit
                printHelpInfo()
                sys.exit()
            else:
                raise AssertionError("unhandled option")

        optList = []
        for opt in opts:
            optList.append(opt[0])

        if (
            ("--qs" in optList)
            and ("--qm" in optList)
            or ("--ts" in optList)
            and ("--tm" in optList)
        ):
            print(
                "\nWarning: the atom number of metal is specified for a metal site. It may be ignored if not coincide with a metal of the site."
            )

        if (not inputType1) or (not inputType2):
            print(
                "\nThe program needs at least two input files to start alignment.\nFor help use -h or --help\n"
            )
            sys.exit()

        else:
            input1 = {"type": inputType1, "path": pathToInput1, "metal": metalID1}
            input2 = {"type": inputType2, "path": pathToInput2, "metal": metalID2}
    else:
        print(
            "\nInput parameters for both input files are mandatory to start alignment.\nFor help use -h or --help\n"
        )
        sys.exit()

    if not args:
        output = None

    elif len(args) > 1:
        usage()
        sys.exit()
    else:
        output = args[0]

    return input1, input2, output, maxDist, rmPdb, scoreOnly

def summary():
    print("\nSummary:\n--------")
    print(
        "MetalS2 is a new algorithm for aligning metal-binding sites based on their three-dimensional structural data.\nThe main application of the tool is detection of structural and functional similarities of sites from protein and nucleotide molecules.\n"
    )


def usage():
    print("\nUsage:\n------")
    print(
        "$./metals2.py [input parameter] <file1> [input parameter] <file2> [input options] <output directory>\n"
    )


def example():
    print("\nExamples:\n------")
    print("$./metals2.py --qs ./data/1ffy_2.site.pdb --ts ./data/1en7_5.site.pdb\n")
    print("or\n")
    print(
        "$./metals2.py --qp ./data/1FFY.pdb --qm 1001 --tp ./data/1EN7.pdb --tm 401\n"
    )


def helpinfo():
    print("Options:\n------")
    print(
        """
     --qp <file1>           input parameter     specify the path to a file with a query pdb
     --tp <file2>           input parameter     specify the path to a file with a target pdb
     --qs <file1>           input parameter     specify the path to a file with a query site
     --ts <file2>           input parameter     specify the path to a file with a target site

     --qm <number>          input option        specify a sequence number of a metal of interest in the query structure
     --tm <number>          input option        specify a sequence number of a metal of interest in the target structure
     --rm_pdb               flag                remove generated PDB and PyMOL files from the output
     --score_only           flag                write only a renamed score file for each alignment

     -d   <number>          input option        specify the maximum distance between atoms to considere two atoms as possible neighbours (in A)

     -h   --help            flag                print help information
     -u   --usage           flag                print usage summary
    """.strip()
    )


def moreinfo():
    print("\nMore info:\n----------")
    print(
        "Input parametes are mandatory.\nInput options, options and output directory are non-mandatory.\nIf output directory is not supplied, results will be stored in a current working directory.\n"
    )


def opened(fileName):
    try:
        f = open(fileName)
        f.close()
        return True
    except OSError:
        return False
