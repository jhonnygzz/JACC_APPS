#pragma once

// Whether or not we retrieved the state of the repo.
#define MINIBUDE_VCS_RETRIEVED_STATE true

// The SHA1 for the HEAD of the repo.
#define MINIBUDE_VCS_HEAD_SHA1 "d19e6eaa36f8a621737655ecb7035ad35ebe8c26"

// Whether or not there were uncommited changes present.
#define MINIBUDE_VCS_IS_DIRTY true

// The name of person who committed HEAD.
#define MINIBUDE_VCS_AUTHOR_NAME "alexishuante"

// The email addess of the person who committed HEAD.
#define MINIBUDE_VCS_AUTHOR_EMAIL "alexis.huante@gmail.com"

// When HEAD was committed.
#define MINIBUDE_VCS_COMMIT_DATE_ISO8601 "2024-10-27 08:20:27 -0400"

// The subject line for the HEAD commit.
#define MINIBUDE_VCS_COMMIT_SUBJECT "Used bm2 (long benchmark)instead of bm1 (short benchmark)"

// The notes attached to the HEAD commit.
#define MINIBUDE_VCS_COMMIT_BODY Small adjustments.

// The output from git --describe (e.g. the most recent tag)
#define MINIBUDE_VCS_DESCRIBE "d19e6ea"
