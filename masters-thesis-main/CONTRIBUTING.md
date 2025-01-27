# Contributing

To contribute to the project, one should follow the guidelines listed in this
file. Otherwise contributions might not be accepted by the project
maintainers.

## Use the ticket system

This project is most likely housed on some Git-related version control
website, which means that there *will* be a ticket system available. If you
notice a problem or come up with an enhancement, please create a ticket in the
ticket system, before doing anything else.

## Create pull/merge requests based on tickets

Both GitHub and GitLab support creating [branches][branch] on issue pages. If
you spot a ticket you would like to work on, mark yourself as the ticket
*assignee* on the issue page, in addition to creating a branch on the same
page. The importance of creating branches on the version control website
instead of on your local copy of the repository is that the website will
automatically give the branch a name that matches the issue or ticket title.
This again makes it easier to keep track of which branch corresponds to which
ticket.

With the branch created, you should create a [pull request] towards the `main`
branch if on GitHub, [merge request] if on GitLab, or something else that
corresponds to these concepts, if on some other website. **Note** that you
might not be able to create a pull request, unless you have already made some
commits on the new branch.

You can keep developing and making commits on the `ticket-branch`, until you
think the `ticket-branch` is ready to be merged into the `main` branch. When
the time comes for a merge…

[branch]: https://git-scm.com/book/en/v2/Git-Branching-Branches-in-a-Nutshell
[pull request]: https://docs.github.com/en/pull-requests
[merge request]: https://docs.gitlab.com/ee/user/project/merge_requests/

## Squash ticket branches before merging

It is not necessarily a good practice to make a million commits on a ticket
branch, and then merge all of those separate commits into the main branch, as
it tends to clutter up project history. It is much cleaner to *squash* the
commits on the ready-to-be-merged `ticket-branch` into a single commit, and
then merge the squashed branch *onto* `main`. The goal here is that the
history of the `main` branch remains linear, which makes it easy to follow.

The squashing is done with an [interactive rebase][rebase]. To put it simply,
run

	git rebase main ticket-branch

so that the commits of the `ticket-branch` are moved linearly after the `main`
branch, resolving any [merge conflicts] that appear along the way. Once the
commits have been moved and conflicts have been resolved, run

	git rebase -i main

This will open up a text editor (see the Git setting `core.editor`
[**here**][git-config]), where you should mark all commits except the topmost
one as `squash`. The topmost row should remain as a `pick` commit, into which
everything else will be squashed:

	pick   commit-hash-1 Commit message 1.
	squash commit-hash-2 Commit message 2.
	squash commit-hash-3 Commit message 3.
	squash commit-hash-4 Commit message 4.
	 ⋯

Once this view is saved and closed, Git will open up a second text editor to
allow you to edit the commit message of the squashed commit. Please retain the
commit messages of the squashed commits during the squash, and add a new title
to the squashed commit, indicating in a shortened manner what changes took
place:

	New title

	Text explaining that a squash happended.

	----------------------------

	Commit message 1.

	Commit message 2.

	Commit message 3.

	Commit message 4.

Once this is done, force push the commit into the remote branch with the  `-f`
switch. The squashed commit is now ready to be merged. Please do this by
closing the pull or merge request on the VCS website.

**Note:** the rebasing (squashing) of branches is considered to be a form of
*rewriting history* in Git. Therefore it should never be performed on the
stable `main` branch, so no established work is accidentally lost. Ticket
branches can be rebased to ones hearts content, as long as care is taken to
make sure no work is lost along the way.

[rebase]: https://www.atlassian.com/git/tutorials/rewriting-history/git-rebase
[merge conflicts]: https://www.atlassian.com/git/tutorials/using-branches/merge-conflicts
[git-config]: https://git-scm.com/book/en/v2/Customizing-Git-Git-Configuration

## Create Git tags for new releases

When a feature is added or a bug is fixed, the commit where this occurred
should be tagged with an *annotated* tag via Git, to signify that the commit
is one of the versions a user should concern themselves with
([documentation][git-tag]). The command to do this with is

	git tag -a vA.B.C

where the `A` is the *major* version, `B` is the *minor* version and `C` is
the *patch* version. The number `A` should be incremented, if the new version
introduces a backwards incompatible change, number `B` should be incremented
if a new backwards-compatible feature was added, and `C` should be
incremented, if a bug in the code was fixed.

When the annotated tagging command is issued, Git will open up a text editor,
where you *must* insert bullet points regarding what changes were made to the
project since the last tagged commit:

	Versio A.B.C

	- Some change.

	- Some other change.

If this tag message is not descriptive enough, the tag *will* be rejected by
the project team.

[git-tag]: https://git-scm.com/book/en/v2/Git-Basics-Tagging
