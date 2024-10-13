import sys

def input_escape(question):
    answer = input(f"{question}\n(press ESCAPE to exit)")
    if answer == "":
        if yes_or_no("EXIT PROGRAM?"):
            sys.exit("PROGRAM TERMINATED BY USER")
    return answer

def select_one(question, options):
    return_key = True
    if not isinstance(options, dict):
        return_key = False
        options = {str(i+1): option for i, option in enumerate(options)}
    pretty_choices = ' -- '.join(': '.join(item) for item in sorted(options.items()))
    choice = input_escape(
        f"""{question}
        {pretty_choices}""")
    while choice not in options:
        choice = input_escape(f"""Invalid choice, possibilities are:\{pretty_choices}\n""")
    return choice if return_key else options[choice]

def yes_or_no(question):
    answer = input_escape(f"""{question} (y/n)\n""").upper()
    while answer not in ["Y", "N"]:
        answer = input_escape(f"""{question} (y/n)\n""").upper()
    return answer == "Y"


def input_list(question):
    return input_escape(f"""{question}\n (input values separated by comma)""").replace(" ", "").split(",")
    