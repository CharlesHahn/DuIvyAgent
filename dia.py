"""
DuIvyAgent

"""


import json

import yaml
import requests

from prompt import system_prompt
from terminal import run_terminal



class NEO(object):

    def __init__(self):
        with open("config.yaml", "r") as fo:
            self.config = yaml.safe_load(fo)
        print(self.config)
        self.messages = [{"role": "system", "content": system_prompt}]
    

    def __call__(self):

        headers = {"Authorization": f"Bearer {self.config['api_key']}", "Content-Type": "application/json"}
        data = {
            "model": self.config['model'],
            "temperature": self.config['temperature'],
            # "max_tokens": self.config['max_tokens'],
            # "top_p": self.config['top_p'],
            # "stream": self.config['stream'],
        }


        user_prompt = input("\n>>>>>>>>>\nUser >>> ")
        while True:
            if user_prompt == "exit":
                break

            self.messages.append({"role": "user", "content": user_prompt})
            data["messages"] = self.messages
            # print(data)
            response = requests.post(
                url=self.config['base_url'],
                headers=headers,
                json=data,
                # timeout=self.config["timeout"]
            )
            response.raise_for_status()  # 自动处理 HTTP 错误
            result = response.json()
            # print(result)
            resp_content = result["choices"][0]["message"]["content"][8:-3]
            self.messages.append(result["choices"][0]["message"])
            # print(self.messages)
            print(f"\nNEO >>> {resp_content[1:-1]}")

            ## parse command
            resp_json = json.loads(resp_content)
            if "command" in resp_json:
                print(f"Are you sure to run >> {resp_json["command"]}", end = " (y/n) ")
                if input().strip() == "y":
                    cmd_res = run_terminal(resp_json["command"])
                    print(f"returncode: {cmd_res["returncode"]}")
                    print(f"output: >>> \n {cmd_res["output"]}")
                    print(f"error: >>> \n {cmd_res["error"]}")
                    user_prompt = str(cmd_res)
                else:
                    user_prompt = input("\n>>>>>>>>>\nUser >>> ")
                    if user_prompt == "exit":
                        break
            else:
                user_prompt = input("\n>>>>>>>>>\nUser >>> ")
                if user_prompt == "exit":
                    break







def main():
    NEO()()



if __name__ == "__main__":
    main()