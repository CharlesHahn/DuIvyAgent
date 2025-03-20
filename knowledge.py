""" 

"""

import json
import os


class KnowledgeBase:
    def __init__(self, file_path: str = "knowledge.json"):
        self.file_path = file_path
        self.knowledge = {}
        self.load_knowledge()
    
    def load_knowledge(self):
        """从JSON文件加载知识库"""
        if os.path.exists(self.file_path):
            with open(self.file_path, 'r', encoding='utf-8') as f:
                pairs = json.load(f)
                for pair in pairs:
                    self.knowledge[pair['query']] = pair['answer']
        else:
            self.knowledge = {}

    def save_knowledge(self):
        """保存知识库到JSON文件"""
        pairs = [{'query': k, 'answer': v} for k, v in self.knowledge.items()]
        with open(self.file_path, 'w', encoding='utf-8') as f:
            json.dump(pairs, f, ensure_ascii=False, indent=2)

    def add_pair(self, query, answer_type, content):
        """添加新条目（自动合并重复query）"""
        # 合并策略：直接覆盖已有条目
        # TODO 合并策略
        self.knowledge[query] = {
            'type': answer_type,
            'content': content
        }
        self.save_knowledge()

    def query_pair(self, query):
        """查询条目（完全匹配）"""
        pair = self.knowledge.get(query)
        if not pair:
            return None
            
        if pair['type'] == 'file':
            # 读取文件内容
            file_path = os.path.dirname(self.file_path) + '/' +  pair['content']
            if os.path.exists(file_path):
                with open(file_path, 'r', encoding='utf-8') as f:
                    return f.read()
            else:
                raise FileNotFoundError(f"知识库记录的文件不存在: {file_path}")
        else:
            # 直接返回字符串内容
            return pair['content']
    
    def get_all_queries(self):
        """返回所有条目的query"""
        return list(self.knowledge.keys())




# 使用示例
if __name__ == "__main__":
    # 初始化知识库（自动创建空json文件）
    kb = KnowledgeBase("knowledge.json")
    
    # 查询条目
    print(kb.query_pair("力场选择"))  # 输出字符串内容
    print(kb.query_pair("md.mdp模板"))  # 输出文件内容
    
    # 更新已有条目
    kb.add_pair(
        query="力场",
        answer_type="string",
        content="推荐使用amber99sb-ildn或charmm36"
    )
    print(kb.query_pair("力场"))  # 显示更新后的内容