# Changelog

This file lists the changes that have been made to this project along its
development. Note that as of version `2.0`, the changes should also be
discoverable via the [Git tags of the project][project-tags-page].

[project-tags-page]: https://gitlab.com/tuni-official/thesis-templates/masters-thesis/-/tags

###### v2.8 Include listings package and uncomment Finnish code example

- Allows the inclusion of code snippets via the command `\lstinputlisting`.
  This command is absolutely required, if the code or its comments contain
  special Unicode characters.

- The package is loaded in `preamble.tex`.

- The Finnish listing example was uncommented in `tex/esitystyyli.tex`, made
  possible by instructing `listings` to convert special characters to
  equivalent LaTeΧ commands in the preamble.

###### v2.7 Uncomment inclusion of siunitx and add \norm command

- There is no point in commenting out the inclusion of the package `siunitx`,
  as it does not increase compilation times noticeably. This also allowed the
  SI-unit example in `tex/esitystyyli.tex` to be uncommented.

- Added the `\norm` command to `preamble.tex`, which allows typesetting the
  norm of a vector 𝐯 as ‖𝐯‖. This is utilized in the now-uncommented SI-unit
  example.

- Added some space before punctuation at the ends of equation examples. This
  makes the punctuation more noticeable.

- Removed an unnecessary hyphen in a paragraph.

###### v2.6.3 Update Usage section in README

- It now contains a link to the canonical Tags-page of the project.

###### v2.6.2 Update Bugs and improvements section in README

- There is now a link to the canonical project Issues page.

- Also added a suggestion to sign up, in case one wants to post issues.

- Since GitLab does not provide a forum for project discussions, instructions
  regarding it were removed.

###### v2.6.1 Change example theorem environment lause → theorem

- Fix the incorrectly set English theorem environment example in the file tex/esitystyyli.tex.

###### v2.6 Define new theorem environments + \code

- The file preamble.tex now defines new theorem environments using the package
  amsthm:

  - maaritelma,
  - lause,
  - apulause,
  - seurauslause and
  - esimerkki.

  These have the English equivalents

  - definition,
  - theorem,
  - lemma,
  - corollary and
  - example.

  These all subscribe to the same counter.

- Added a (Finnish) section under the writing style chapter, on how these
  environments should be utilized.

- Also defined the command \code for inline teletyped text.

###### v2.5.3 Fix load order of pdfx + improve main.tex

- The loading of the package pdfx is moved to the end of tauthesis.cls, in
  order to facilitate the setting of the title language in the metadata file
  \jobname.xmpdata, before the information is written to the output PDF.

- Move PDF version setting to the file set-pdf-version.tex, to again make the
  file main.tex more readable. This is not something that the user should
  worry about, so this should again provide for a better user experience.

###### v2.5.2 Update tauthesis.cls date and CHANGELOG

- Updates the date specified in tauthesis.cls to be in line with the tags. The
  file version now follows the project tags.

- Brings the file CHANGELOG.md up to date with the recently added Git-tags.

- Also unifies the format between them, so that things are easy to copy over.

###### v2.5.1 Move chapter and label commands to respective chapter files

- Move \chapter and the respective \label commands into the files that contain
  the chapters, as it makes more sense for everything related to a chapter to
  be in the file that contains it.

- Added comment symbols % at the ends of the \chapter command lines, before
  the \label commands, so that any TeΧ compilers would interpret the \label
  commands as being on the same line with each \chapter. This might get rid of
  some warnings related to page numbering.

###### v2.5 Move preamble and title page to their own files

- Moved preamble from main.tex → preamble.tex.

- Moved title page from main.tex → titlepage.tex.

- Automated primary language selection on title page via language commands
  defined in the file tauthesis.cls. Now a user no longer needs to modify
  these themselves, which should provide for a better user experience.

- Added metadata commands \myyear, \mymonth and \myday in main.tex.

- Changed the command invocation \documentclass{tauthesis} →
  \documentclass[finnish]{tauthesis}, to make it clear how the language is
  selected.

###### v2.4.2 Update © and maintainer info

- Copyright year from 2018 → 2023.

- Maintainer from Ville Koljonen → Santtu Söderholm.

###### v2.4.1 Update documentation

- Added short-circuiting to the compilation sequence in the README, so that
  later commands in the sequence are not attempted if a previous one fails
  with an exit code other than 0.

- Also marked the related code block as a shell script, so that it would be
  highlighted when the README is viewed with a suitable reader.

###### v2.4 Replace amsmath with mathtools

- Replaces the included package amsmath with its patched and extended version,
  called mathtools.

###### v2.3.2 Utilize subtitle commands from metadata section

- Fixed the subtitle defined by the user in the metadata section of main.tex
  not being displayed on the title page of the thesis.

###### v2.3.1 Add instructions related ot Git-tagging

- Add instructions related to Git tagging to CONTRIBUTING.md.

###### v2.3 Simplify metadata insertion

- Simplify metadata insertion by defining metadata via commands that get
  reused.

- This metadata section is close to the start of main.tex, and users should
  simply replace each myvalue in the commands \def\mykey{myvalue}.

###### v2.2 Allow building the project with LuaLaTeX

- Allows building the project with lualatex.

- Note that the axessibility package is unavailable, if lualatex is used.

###### v2.1 Update formatting in Finnish lists of references + CHANGELOG and CONTRIBUTING files

- Change "et al." → "ym." in Finnish lists of references.

- Add separate CHANGELOG.md and CONTRIBUTING.md files.

- Update README with more detailed usage instructions.

###### v2.0.1 Add a .gitignore file

- Added a .gitignore file to prevent unwanted files from being uploaded into the repository.

###### v2.0 Original upload by Ville Koljonen

- Contains the repository as originally uploaded by Ville Koljonen.
