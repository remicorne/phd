import sys

def input_escape(f):
    def wrapper(*args, **kwargs):
        answer = input(f"{args[0]}\n(press ESCAPE to exit)")
        if answer == "":
            if Question.yes_or_no("EXIT PROGRAM?"):
                sys.exit("PROGRAM TERMINATED BY USER")
            else: return f(*args, **kwargs)
        return answer
    return wrapper
    
class Question:
    
    def input_escape(question):
        answer = input(f"{question}\n(press ESCAPE to exit)")
        if answer == "":
            if Question.yes_or_no("EXIT PROGRAM?"):
                sys.exit("PROGRAM TERMINATED BY USER")
        return answer
    
    def input(question):
        return Question.input_escape(question)
    
    def select_one(question, options):
        if not isinstance(options, dict):
            options = {str(i+1): option for i, option in enumerate(options)}
        choice = Question.input_escape(
            f"""{question}
            {options}""")
        while choice not in options:
            choice = Question.input_escape(f"""Invalid choice, possibilities are:\{options}\n""")
        return options[choice]

    def yes_or_no(question):
        answer = Question.input_escape(f"""{question} (y/n)\n""").upper()
        while answer not in ["Y", "N"]:
            answer = Question.input_escape(f"""{question} (y/n)\n""").upper()
        return answer == "Y"