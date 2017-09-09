#!/bin/sh

cat << EOF | sqlite3 enron.db

.read "schema.sql"
.separator "\t"
.import "employees.tsv" EmployeeBase
.import "messages.tsv" MessageBase
.import "recipients.tsv" RecipientBase

EOF
