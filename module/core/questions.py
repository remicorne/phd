import sys

def input_escape(f):
    def wrapper(*args, **kwargs):
        answer = input(f"{args[0]}\n(press ESCAPE to exit)")
        if answer == "":
            if Questions.yes_or_no("EXIT PROGRAM?"):
                sys.exit("PROGRAM TERMINATED BY USER")
            else: return f(*args, **kwargs)
        return answer
    return wrapper
    
class Questions:
    
    def input_escape(question):
        answer = input(f"{question}\n(press ESCAPE to exit)")
        if answer == "":
            if Questions.yes_or_no("EXIT PROGRAM?"):
                sys.exit("PROGRAM TERMINATED BY USER")
        return answer
    
    def input(question):
        return Questions.input_escape(question)
    
    def select_one(question, options):
        choice = Questions.input_escape(
            f"""{question}
            {options}""")
        while choice not in options.keys():
            choice = Questions.input_escape(f"""Invalid choice, possibilities are:\{options}\n""")
        return choice

    def yes_or_no(question):
        answer = Questions.input_escape(f"""{question} (y/n)\n""").upper()
        while answer not in ["Y", "N"]:
            answer = Questions.input_escape(f"""{question} (y/n)\n""").upper()
        return answer == "Y"