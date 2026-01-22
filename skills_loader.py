"""
Skills Loader Module

This module provides functionality to load and manage AI agent skills from the skills directory.
Skills are organized as folders containing SKILL.md files with YAML frontmatter.
"""

import os
import yaml
from pathlib import Path
from typing import Dict, List, Optional


class SkillsLoader:
    """Load and manage AI agent skills from the skills directory."""

    def __init__(self, skills_dir: str = "skills"):
        """
        Initialize the skills loader.

        Args:
            skills_dir: Path to the skills directory (default: "skills")
        """
        self.skills_dir = Path(skills_dir)
        self.skills: Dict[str, dict] = {}
        self.load_skills()

    def load_skills(self) -> None:
        """Scan the skills directory and load all skills."""
        if not self.skills_dir.exists():
            print(f"Warning: Skills directory '{self.skills_dir}' does not exist.")
            return

        # Recursively find all SKILL.md files
        skill_files = list(self.skills_dir.rglob("SKILL.md"))

        for skill_file in skill_files:
            try:
                skill_info = self._parse_skill_file(skill_file)
                if skill_info:
                    skill_name = skill_info['name']
                    self.skills[skill_name] = skill_info
            except Exception as e:
                print(f"Warning: Failed to load skill from {skill_file}: {e}")

    def _parse_skill_file(self, skill_file: Path) -> Optional[dict]:
        """
        Parse a SKILL.md file and extract metadata and content.

        Args:
            skill_file: Path to the SKILL.md file

        Returns:
            Dictionary containing skill information or None if parsing fails
        """
        with open(skill_file, 'r', encoding='utf-8') as f:
            content = f.read()

        # Extract YAML frontmatter
        frontmatter = self._extract_frontmatter(content)

        if not frontmatter or 'name' not in frontmatter:
            print(f"Warning: Invalid skill file {skill_file}: missing name in frontmatter")
            return None

        # Remove frontmatter from content
        markdown_content = self._remove_frontmatter(content)

        return {
            'name': frontmatter.get('name'),
            'description': frontmatter.get('description', ''),
            'allowed_tools': frontmatter.get('allowed-tools', []),
            'path': str(skill_file),
            'content': markdown_content,
            'full_content': content
        }

    def _extract_frontmatter(self, content: str) -> Optional[dict]:
        """
        Extract YAML frontmatter from markdown content.

        Args:
            content: Full content of the SKILL.md file

        Returns:
            Parsed YAML frontmatter as dictionary or None if not found
        """
        lines = content.split('\n')

        if not lines or lines[0] != '---':
            return None

        # Find the end of frontmatter
        end_idx = -1
        for i in range(1, len(lines)):
            if lines[i] == '---':
                end_idx = i
                break

        if end_idx == -1:
            return None

        # Extract and parse YAML
        frontmatter_str = '\n'.join(lines[1:end_idx])
        try:
            return yaml.safe_load(frontmatter_str)
        except yaml.YAMLError as e:
            print(f"Warning: Failed to parse YAML frontmatter: {e}")
            return None

    def _remove_frontmatter(self, content: str) -> str:
        """
        Remove YAML frontmatter from markdown content.

        Args:
            content: Full content of the SKILL.md file

        Returns:
            Content without frontmatter
        """
        lines = content.split('\n')

        if not lines or lines[0] != '---':
            return content

        # Find the end of frontmatter
        end_idx = -1
        for i in range(1, len(lines)):
            if lines[i] == '---':
                end_idx = i
                break

        if end_idx == -1:
            return content

        # Return content without frontmatter
        return '\n'.join(lines[end_idx + 1:])

    def get_skill_list(self) -> List[str]:
        """
        Get list of available skill names.

        Returns:
            List of skill names
        """
        return list(self.skills.keys())

    def get_skill_metadata(self, skill_name: str) -> Optional[dict]:
        """
        Get metadata for a specific skill.

        Args:
            skill_name: Name of the skill

        Returns:
            Dictionary containing skill metadata or None if not found
        """
        if skill_name not in self.skills:
            return None

        skill = self.skills[skill_name]
        return {
            'name': skill['name'],
            'description': skill['description'],
            'allowed_tools': skill['allowed_tools']
        }

    def get_skill_content(self, skill_name: str) -> Optional[str]:
        """
        Get the full content of a specific skill.

        Args:
            skill_name: Name of the skill

        Returns:
            Full content of the skill (including frontmatter) or None if not found
        """
        if skill_name not in self.skills:
            return None

        return self.skills[skill_name]['full_content']

    def get_skill_markdown(self, skill_name: str) -> Optional[str]:
        """
        Get the markdown content of a specific skill (without frontmatter).

        Args:
            skill_name: Name of the skill

        Returns:
            Markdown content of the skill (without frontmatter) or None if not found
        """
        if skill_name not in self.skills:
            return None

        return self.skills[skill_name]['content']

    def get_all_skills_metadata(self) -> List[dict]:
        """
        Get metadata for all available skills.

        Returns:
            List of dictionaries containing skill metadata
        """
        return [
            self.get_skill_metadata(skill_name)
            for skill_name in self.get_skill_list()
        ]

    def format_skills_for_prompt(self) -> str:
        """
        Format skills information for inclusion in system prompt.

        Returns:
            Formatted string listing available skills
        """
        if not self.skills:
            return "No skills available."

        lines = ["## Available Skills\n"]
        for skill_name in self.get_skill_list():
            skill = self.skills[skill_name]
            lines.append(f"- **{skill_name}**: {skill['description']}")
            lines.append(f"  → Use `get_skill_content('{skill_name}')` to load full instructions")

        return '\n'.join(lines)


# Example usage
if __name__ == "__main__":
    loader = SkillsLoader()

    print("=" * 80)
    print("Available Skills:")
    print("=" * 80)
    for skill_name in loader.get_skill_list():
        metadata = loader.get_skill_metadata(skill_name)
        print(f"\n[{skill_name}]")
        print(f"Description: {metadata['description']}")
        if metadata['allowed_tools']:
            print(f"Allowed tools: {', '.join(metadata['allowed_tools'])}")

    print("\n" + "=" * 80)
    print("Formatted for Prompt:")
    print("=" * 80)
    print(loader.format_skills_for_prompt())
