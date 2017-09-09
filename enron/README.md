The Enron Email Corpus
======================

The Enron email corpus is coded into three tables. The `Employee` table
contains the information about the 156 employees in the dataset. The
`Message` and `Recipient` tables together contain information about the
messages. It is necessary to have two tables since some messages have multiple
recipients.


The next example illustrates how the messages get coded in the tables.
Message 4 was sent by employee 138 to employees 109, 49, 120, and 59. This
gets coded as a single row in the `Message` table and four rows in the
`Recipient` table. The row in the `Message` table is

    mid  filename          unix_time  subject                 sender_id
    ---  ----------------  ---------  ----------------------  ---------
    4    taylor-m/sent/23  911874180  Re: Coral Energy, L.P.  138
    
while the rows in the `Recipient` tables are

    mid  rno  receiver_id
    ---  ---  -----------
    4    1    109
    4    2    49
    4    3    120
    4    4    59


The end of the document contains detailed information about the individual 
tables.


Zhou et al. did the majority of the work in preprocessing the Enron corpus
[1]. Perry & Wolfe extracted the employee genders, seniorities, and
departments from Zhou et al.'s employee information table and the raw message
contents [2]. If you use the data in a publication, please cite the references
below.

  1. Perry, P. O., and Wolfe, P.J. (2011).  "Point process modeling for
     directed interaction networks." _Submitted_.

  2. Zhou, Y., Goldberg, M., Magdon-Ismail, M., and Wallace, W. A. (2007).
     "Strategies for cleaning organizational emails with an application to 
     Enron email dataset." _5th Conf. of North American Association for 
     Computational Social and Organizational Science (NAACSOS 07)_.


SQLite
------

The script `create.sh`, when run from the current directory, creates an
SQlite3 database called `enron.db`. It then imports the data from
`employees.tsv`, `messages.tsv`, and `recipients.tsv`. For the script to work,
the `sqlite3` program must be in your path. See `schema.sql` for the database
schema.

The raw data files from Zhou et al. are `employeeDirsChp3New.txt` and
`EmpInfo.pdf`.


Employee Information (employees.tsv)
------------------------------------

The table of 156 employees. Columns are as follows:

  - employee id (1-156)
  - full name
  - department (Legal, Trading, or Other)
  - long department
  - title
  - gender (Female or Male)
  - seniority (Junior or Senior)
  

Message Information (messages.tsv)
----------------------------------

The table of 21,635 messages. Columns are as follows:

  - message id (1-21635)
  - filename
  - unix time (seconds since January 1, 1970)
  - message subject
  - sender's employee id (1-156)
  

Recipient Information (recipients.tsv)
--------------------------------------

The table of 38,388 recipients for each of the messages. Columns as follows:

  - message id (1-21635)
  - recipient number (1-57)
  - receiver's employee id (1-156)
