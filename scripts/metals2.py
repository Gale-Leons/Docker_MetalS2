#!/usr/bin/env python

import os
import sys

from lib.alignment import alignSites
from lib.reporting import reportAlignment, storeAlignment
from lib.utilites import processCommandLineArguments

if __name__ == "__main__":
    try:
        print("\nPre-processing the input data. Please, wait...")
        # informazioni provenienti dalla riga di comando per il lancio del programma da riga di comando
        # 1 contro tutti
        # maxDist non veniva esplicitato quindi verra' preso il Default
        sitesList1, sitesList2, pathToOutput, maxDist = processCommandLineArguments(
            sys.argv
        )
        print("Ready to go.\n")

        multi = False
        if (
            len(sitesList1) > 1 and len(sitesList2) > 1
        ):  # caso multiplo allineamento > incrociato ?
            if not os.path.exists(pathToOutput):
                os.mkdir(pathToOutput)
            multi = True

        root = pathToOutput
        for site1 in sitesList1:
            if multi:
                pathToOutput = os.path.join(pathToOutput, site1.name)

            for site2 in sitesList2:
                print(f"Aligning {site1.name} with {site2.name}...")
                try:
                    p = "pairwise"
                    d = "database"
                    alignedSite, scoreReport = alignSites(
                        site1, site2, maxDist, p
                    )  # Funzione che svolge l-allineamento

                except Exception as e:
                    print(e)
                    print("The alignment process can't be completed.")
                    continue
                else:
                    sitesDirRoot = storeAlignment(
                        site1, alignedSite, pathToOutput
                    )  # altra funzione
                    reportAlignment(
                        site1, alignedSite, scoreReport, sitesDirRoot, maxDist
                    )  # altra funzione di Report
                    print("Done.\n")

            pathToOutput = root

        print(f"Results are stored here: {pathToOutput}\n")

    except SystemExit:
        pass
