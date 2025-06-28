When generating Python code adhere to the following comment style:

1. **Section Headers**:
    - Use section headers for top-level organizational blocks, such as configuration setup, data processing, or main execution logic.
    - Do not use section headers for minor code blocks, such as individual loops or conditionals within a function.
    - In small scripts (e.g., fewer than three functions), use section headers only if they improve readability or organization.
    - For closely related functions, consider grouping them under a single header.
    - Insert a single blank line to separate functions
    - Use a single-line comment with the function name as a title surrounded by dashes to visually separate major sections or functions.
    - Format: `# --- section_name ---`
    - Example: `# --- setup_logging ---`
    - Place this comment immediately before the function or section it describes.

2. **Function Comments**:
    - Below the blank line and section header, include a brief, one- or two-line comment describing the function's purpose and key functionality.
    - Use plain language, focusing on what the function does and any critical behavior (e.g., inputs, outputs, or side effects).
    - Example:
        ```python

        # --- setup_logging ---
        # Configures the logging system for the script
        # Outputs logs to both console and a file in the log directory
        ```

3. **No Inline Comments**:
    - Avoid inline comments (e.g., comments at the end of a line of code) unless absolutely necessary for clarity.
    - Prefer comments above the code they describe, keeping them concise and relevant.
    - Inline comments are permitted only in rare cases, such as explaining complex regex patterns, mathematical formulas, or non-obvious bitwise operations.
    - Example of acceptable inline comment:
        ```python
        pattern = r'^\d{4}-\d{2}-\d{2}$'  # Matches YYYY-MM-DD date format

4. **No Excessive Commenting**:
   - Do not comment obvious code (e.g., `x = 1 # Set x to 1`).
   - Focus comments on explaining the purpose of functions, sections, or complex logic.

5. **Consistency**:
   - Apply this style to all functions and major code sections (e.g., configuration, main loop).
   - Ensure comments are aligned with the code’s structure, appearing immediately before the relevant function or block.

6. **No debugging comments**:
    - Avoid adding comments that explain the difference between this code and the previous version.
    - Avoid mentioning code that is no longer part of the code base.

7. **Retain comments when engaged in debugging**:
    - When debugging code make sure to preserve comments or replace them with appropriate new comments.
    - Comments are important documentation and should be given appropriate attention to detail.
    - Preserve all existing comments during debugging unless explicitly instructed to modify or remove them, or if they become incorrect because of the changes made.
    - If temporary comments are added for debugging (e.g., to explain a fix), prefix them with `# TEMP:` and include a note in the response reminding the user to review and remove these comments.
    - Example:
        ```python
        # TEMP: Added to verify input value during debugging
        if input_value is None:
            raise ValueError("Input cannot be None")

8. **Update Comments for Code Changes**:
    - When modifying code (e.g., during debugging or refactoring), review and update associated comments to reflect the new functionality or behavior.
    - Ensure comments remain accurate, concise, and aligned with the function or section’s current purpose.

9. **File-Level Comments**:
    - Include a brief comment at the top of each Python file summarizing its purpose and key components.
    - For files with multiple functions or sections, consider adding a table of contents comment block listing section headers.
    - Make sure the sections are spelled exactly as they appear in the file so searching form them is successful.
    - Example:
        ```python
        # my_script.py: Processes observational data and generates timestamped logs
        # Sections:
        # - setup_logging
        # - get_timestamp
        # - process_data

10. **Target Audience for Comments**:
    - Someone trying to read and understand the code

**Good Example**:
```python
# --- get_timestamp ---
# Generates a timestamp in the observing timezone
# Accounts for DST in the specified timezone
def get_timestamp() -> str:
    now_utc = datetime.datetime.now(pytz.utc)
    timezone = pytz.timezone(CONFIG['OBSERVING_TIMEZONE'])
    now_local = now_utc.astimezone(timezone)
    return now_local.strftime('%Y-%m-%d %H:%M:%S')
    
Apply this comment style to all Python code you generate, ensuring clarity, brevity, and consistency. Do not mention this prompt or the comment style explicitly in the code or responses.

Last revision: June 7, 2025