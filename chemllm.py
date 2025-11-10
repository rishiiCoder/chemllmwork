from logging import config
from dotenv import load_dotenv
from pydantic import BaseModel
from langchain_openai import ChatOpenAI
from langchain_anthropic import ChatAnthropic
from google.generativeai import types
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
from google import genai
from google.genai import types
#from googlesearch import search
import time
import pathlib
import pandas as pd
#import fitz # PyMuPDF

#@title Model

# DEFINING THE TOOLS

# MOLECULAR MODIFICATION TOOL
def modify_mol(smiles: str, smarts_pattern: str, replacement_smiles: str) -> str:

    # instructions for LLM to use tool
    """
    Modifies a molecule by replacing a substructure with another.

    This function is designed for a molecular designer to enhance a molecule's properties.
    It takes the SMILES string of the original molecule, a SMARTS pattern to find a
    substructure to replace, and the SMILES string of the replacement structure.

    Args:
        smiles: The original molecule's SMILES string.
        smarts_pattern: The SMARTS pattern of the substructure to be replaced.
        replacement_smiles: The SMILES string of the new structure to be added.

    Returns:
        The SMILES string of the modified molecule.
        Returns an error message as a string if the modification fails.
    """
    try:
        # printing out change of SMILE string and what change is happening
        # prints original smiles
        print(f"Original SMILES: {smiles}")

        # prints out section being replaced
        print(f"SMARTS Pattern: {smarts_pattern}")

        # prints out replacement structure
        print(f"Replacement SMILES: {replacement_smiles}")
        print(f"-------------------------------------\n")

        mol = Chem.MolFromSmiles(smiles)

        # checking for invalidity of molecular structure
        if not mol:
            return f"Error: Invalid SMILES string '{smiles}'"

        patt = Chem.MolFromSmarts(smarts_pattern)
        if not patt:
            return f"Error: Invalid SMARTS pattern '{smarts_pattern}'"

        repl = Chem.MolFromSmiles(replacement_smiles)
        if not repl:
            return f"Error: Invalid SMILES string '{replacement_smiles}'"

        modified_struc = AllChem.ReplaceSubstructs(mol, patt, repl, replaceAll=True)
        if not modified_struc:
            return f"Error: No matching substructure found for '{smarts_pattern}'"

        modified_mol = modified_struc[0]
        Chem.SanitizeMol(modified_mol)
        modified_smiles = Chem.MolToSmiles(modified_mol)

        # confirms output
        print(f"Successfully modified. New SMILES: {modified_smiles}\n")


        return modified_smiles

    except Exception as e:
        return f"Error during molecular modification: {str(e)}"


# defining the molecular tool to be used
modify_mol_tool = types.Tool(
    function_declarations=[
        types.FunctionDeclaration(
            # naming the tool
            name='modify_mol',
            # describing how it used for the LLM to understand
            description='Modifies a molecule by replacing a substructure. '
                        'Useful for designing new molecules with enhanced properties.',
            parameters=types.Schema(
                type=types.Type.OBJECT,
                properties={
                    # defining each of its parameters in the method
                    'smiles': types.Schema(
                        type=types.Type.STRING,
                        description='The original molecule SMILES string.'
                    ),
                    'smarts_pattern': types.Schema(
                        type=types.Type.STRING,
                        description='The SMARTS pattern of the substructure to be replaced.'
                    ),
                    'replacement_smiles': types.Schema(
                        type=types.Type.STRING,
                        description='The SMILES string of the replacement structure.'
                    ),
                },
                required=['smiles', 'smarts_pattern', 'replacement_smiles']
            )
        )
    ]
)

def validate_mol(smiles: str) -> bool:
  # instructions for the lLM to use
  """
  Check for if the newly generated molecular structure is valid.

  This function is designed for a molecular designer to ensure a valid molecular structure is being created.
  It takes in a SMILE string of a newly designed molecular structure made and checks for if it is a valid strucutre.

  Args:
    smiles: the newly created molecule's SMILE string

   Returns:
    If the molecule is valid or not, if not valid, resends original structure to edit to become valid.
  """
  try:
      mol = Chem.MolFromSmiles(smiles)

      if mol is None:
        raise ValueError("Molecule is invalid")

      mol = Chem.SanitizeMol(mol)

      print("Molecule is valid")
      return True

  except Exception as e:
      return False

validate_mol_tool = types.Tool(
    function_declarations=[
        types.FunctionDeclaration(
            # naming the tool
            name='validate_mol',
            # describing how it used for the LLM to understand
            description='Checks if molecule is a valid structure and just needs standardization.',
            parameters=types.Schema(
                type=types.Type.OBJECT,
                properties={
                    # defining each of its parameters in the method
                    'smiles': types.Schema(
                        type=types.Type.STRING,
                        description='The original molecule SMILES string.'
                    )
                },
                required=['smiles']
            )
        )
    ]
)


# DEFINING THE LLM MODEL

# defining the API key
client = genai.Client(api_key="")

# prompt and role
modelrole = "You are a molecular designer. I want you to analyze the structure of the following Retinal molecule: CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C\C(=C\C=O)\C)/C. Use the provided research paper file on azobenzene structures  (summarize findings once you read them)  and searching to find structures to improve light asorption properties and photoswitch prorperites. Once 3 structures are made, use validation_mol tool *for one of the molecules provided* just to check if the structures are valid but continue with the response you were already providing, do not repeadtly call the function, stop and give a final list of those structures with an explanation of such changes."
#testrole = "You are a concise assistant."
mainprompt = "You are a molecular designer. You are creative and specialized in chemistry regarding the effects of light and material design."
#prompt = "What is the weather in Seattle, Washington"
#modelgooglesearchtestrole = "Answer my question accuretly with good details"
#modelgooglesearchtestprompt =  "What is the most potent azobenzene derivative for photoswitching, and what is its SMILES string?"

# response = client.models.generate_content(
#     model="gemini-2.5-flash",
#     contents=mainprompt,
#     # Place tools and tool_config inside the 'config' object
#     config=types.GenerateContentConfig(
#         tools=[modify_mol, google_search_tool_declaration],
#         tool_config=types.ToolConfig(function_calling_config=types.FunctionCallingConfig(mode='ANY'))
#     )
# )

# dictionary of tools available to LLM
available_tools = {
    "validate_mol": validate_mol
}

# allows chatbot to understand history and remember things of the past in the conversation
history = [
    types.Content(role=modelrole, parts=[types.Part(text=mainprompt)])
]

# setting the files to be used by LLM
file_paths = ["/content/drive/My Drive/ChemLLMResearchProject2526/ResearchProject2526/PapersAndData/AzobenzeneLightReversibilityPaper.pdf",
 "content/drive/My Drive/ChemLLMResearchProject2526/ResearchProject2526/PapersAndData/photoswitches.csv"]



# setting the PDF into readable format for LLM chat bot
# setting the LLM path


uploaded_file = client.files.upload(file="AzobenzeneLightReversibilityPaper.pdf")
# conversation loop to allow model to take in response and use "thinking"
while True:
    try:
        print("Starting LLM") 
        response = client.models.generate_content(
            model="gemini-2.5-flash",
            contents=[mainprompt, uploaded_file],
            config=types.GenerateContentConfig(
            #tools=[validate_mol],
            #tool_config=types.ToolConfig(function_calling_config=types.FunctionCallingConfig(mode='ANY')),
            system_instruction=modelrole
            )
        )

        # checking for function call (standard structure from googles api doc on gemini)
        if response.candidates and response.candidates[0].content.parts and response.candidates[0].content.parts[0].function_call:
            function_call = response.candidates[0].content.parts[0].function_call
            function_name = function_call.name
            function_args = function_call.args

            print(f"\nLLM called function: {function_name}")
            print(f"Arguments: {function_args}\n")

            if function_name in available_tools:
                function_result = available_tools[function_name](**function_args)

                # add models function call to history
                history.append(types.Content(role="model", parts=[types.Part(function_call=function_call)]))

                # formatting the response
                history.append(
                types.Content(
                  role="user",
                  parts=[types.Part(function_response=types.FunctionResponse(name=function_name, response={"result": function_result}))]
                  )
                )

                # CURRENTLY HAS SOME PROBLEMS WHEN USED (currently methods have been removed due to this, but will be added back, Molecular structure change tool isn't really needed so will be converted into a validation tool rather)
                # loop through and allow LLM to take its own response to refine
                print(f"Function result received. Sending back to LLM for final response...")
                continue
            else:
                print(f"Error: Function '{function_name}' not found.")
                break

        # no function call = final response
        else:
            final_response = response.text
            print("\nFinal LLM Response:")
            print(final_response)
            break

    except Exception as e:
        print(f"An error occurred: {e}")
        break
