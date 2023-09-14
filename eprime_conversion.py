import re
import os
import pandas as pd
import numpy as np
import math


class eprime_data:
    """read_eprime.
    convert into a `pandas.DataFrame` behavioral data from eprime txt file.

    Parameter
    ---------
    filename : (`str`)
    The eprime txt file name.

    Attribute
    ---------
    dataframe : (`pandas.DataFrame`)

    """

    def __init__(self, filename):
        self.filename = filename

        # Read txt file
        with open(self.filename, "r", encoding="UTF-16") as f:
            lines = f.readlines()

        # Prepare a list for all entries in the txt file
        rawentries = list()

        # Append to the list the entries read in the file

        for entries in lines:
            e = [re.sub(exp, "", entries) for exp in [r"\s", r"\n", " "]][0]
            e.replace(" ", "")
            rawentries.append(e)

        # find all occurence (indices) of the level 3 log (trials)
        indices_level3 = [
            i for i in range(len(rawentries)) if "Level:3" in rawentries[i]
        ]

        # Prepare a dataframe to store the data
        df = pd.DataFrame(
            columns=[
                "Block",
                "PreTargetDuration",
                "TrialNumber",
                "WarningType",
                "CueType",
                "FlankerType",
                "TargetPosition",
                "TargetImage",
                "Response",
                "CorrectAnswer",
                "IntervalBetweenCueAndTarget",
                "TargetType",
                "TargetDirection",
                "FixationStartOnsetTime",
                "WarningOnsetTime",
                "TargetOnsetTime",
                "ReactionTime",
            ]
        )

        # Loop over all trials (Level 3)

        for i in range(len(indices_level3)):
            # Find the start and end of the trial
            startTrialIdx = indices_level3[i] + rawentries[indices_level3[i] :].index(
                "***LogFrameStart***"
            )
            stopTrialIdx = indices_level3[i] + rawentries[indices_level3[i] :].index(
                "***LogFrameEnd***"
            )

            # Extract the entries for the trial
            pertrial = rawentries[startTrialIdx + 1 : stopTrialIdx]
            trialentries = {
                str.split(entry, ":")[0]: str.split(entry, ":")[1] for entry in pertrial
            }

            # Extract the data when trial "Procedure"  are not practice

            if trialentries.get("Procedure") == "TrialProc":
                if trialentries["TargetDirection"] == "left":
                    direction = 1
                elif trialentries["TargetDirection"] == "right":
                    direction = 2

                # Define correct, incorrect and not answered responses

                if (
                    trialentries["SlideTarget.RESP"] != ""
                    or 0 < int(trialentries["SlideTarget.RT"]) < 1700
                ):
                    if int(trialentries["SlideTarget.RESP"]) == direction:
                        CorrectAnswer = "correct"
                    elif int(trialentries["SlideTarget.RESP"]) != direction:
                        CorrectAnswer = "incorrect"
                elif (
                    trialentries["SlideTarget.RESP"] == ""
                    or trialentries["SlideTarget.RESP"] == "not_answered"
                    or int(trialentries["SlideTarget.RT"]) >= 1700
                    or int(trialentries["SlideTarget.RT"]) == 0
                ):
                    CorrectAnswer = "not_answered"
                    trialentries["SlideTarget.RT"] = 0

                # Define duration of fixation
                FixationStartTime = int(trialentries["SlideFixationStart.OnsetTime"])
                TargetOnsetTime = int(trialentries["SlideTarget.OnsetTime"])
                PreTargetDuration = TargetOnsetTime - FixationStartTime

                # Put in a dictionary the data during the trial
                dfentries = {
                    "Block": int(trialentries["TrialList.Cycle"]),
                    "TrialNumber": int(trialentries["TrialList.Sample"]),
                    "WarningType": trialentries["WarningType"],
                    "CueType": "spatial"

                    if trialentries["WarningType"] == "down"
                    or trialentries["WarningType"] == "up"
                    else trialentries["WarningType"],
                    "FlankerType": trialentries["FlankerType"],
                    "TargetImage": trialentries["TargetImage"],
                    "Response": trialentries["SlideTarget.RESP"],
                    "CorrectAnswer": CorrectAnswer,
                    "IntervalBetweenCueAndTarget": int(
                        trialentries["IntervalBetweenCueAndTarget"]
                    ),
                    "PreTargetDuration": PreTargetDuration,
                    "TargetType": "dwn"

                    if trialentries["TargetType"] == "down"
                    else trialentries["TargetType"],
                    "TargetDirection": trialentries["TargetDirection"],
                    "FixationStartOnsetTime": int(
                        trialentries["SlideFixationStart.OnsetTime"]
                    ),
                    "ReactionTime": int(trialentries["SlideTarget.RT"]) if int(trialentries["SlideTarget.RT"]) != 0 else np.nan,
                }
                # Append the dictionary to the dataframe
                df.loc[i] = dfentries

        # Set the index of the dataframe to the trial number
        df = df.set_index(["TrialNumber"])
        df["CueType"] = df["CueType"].replace("no", "no_cue")
        self.dataframe = df

    def export2csv(self, out_filename=None, **kwargs):
        directory = os.path.dirname(out_filename)

        if not os.path.isdir(directory):
            os.mkdir(directory)
        self.dataframe.to_csv(path_or_buf=out_filename, **kwargs)
