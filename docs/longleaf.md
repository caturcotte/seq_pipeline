## Getting an account on Longleaf
From https://help.rc.unc.edu/request-a-cluster-account:
Faculty, staff, graduate and undergraduate students can request an account. If you have a faculty PI, adviser or lab head, please note that in the sponsor field or comments.

Request an account by following these steps:

1. Go to help.unc.edu and log in with your onyen.
2. Search for Compute Cluster Access in the catalog section of help.unc.edu.
3. Click on Computer Cluster Access.
4. Complete the service request form and click “Request”.

You will be notified by email when your account is ready.

**NOTE:** If the search of help.unc.edu for the Compute Cluster Access service request turns up empty and you are sure you are logged into help.unc.edu with your onyen, please send email to research@unc.edu with a screen shot of your entire help.unc.edu search results page, including the top banner portion and indicate that you cannot request the account you need. Answer the following in that email:
- Your name and onyen:
- Your Department:
- Phone number to reach you while your jobs are running (if we need to):
- Your Email address:
- Your preferred shell (bash or tcsh):
- Faculty sponsor name:
- Faculty sponsor onyen:
- Type of subscription: Longleaf or dogwood
- A description of the type of work you will do on the cluster

## Logging in

[[x11]] forwarding (required for graphical displays such as [[Integrated Genome Viewer|IGV]]:
```
ssh -X onyen@longleaf.unc.edu
```

This should work automatically on Linux. For other systems:

**Windows:**

Windows users should download MobaXterm (Home Edition). Then use the Session icon to create a Longleaf SSH session using longleaf.unc.edu for “Remote host” and your ONYEN for the “username” (Port should be left at 22).

**Mac:**

Mac users can use ssh from within their Terminal application to connect to Longleaf. Be sure to use your UNC ONYEN and password for the login:

```
ssh -X onyen@longleaf.unc.edu
```

To enable x11 forwarding Mac users will need to download, install, and run Xquartz on their local machine in addition to using the “–X” ssh option. Furthermore, in many instances for x11 forwarding to work properly Mac users need to use the Terminal application that comes with Xquartz instead of the default Mac terminal application.

## Main Directory Spaces

NAS home space

Your home directory will be in /nas/longleaf/home/onyen and is backed up via snapshots.

Your home directory has a quota which you will want to monitor occasionally: 50 GB soft limit and a 75 GB hard limit.

### /work storage

Your /work directory will be: /work/users/o/n/onyen (the “o/n/” are the first two letters of your ONYEN) with a quota of 30 TB.

/work is built for high-throughput and data-intensive computing, and intended for data actively being computed on, accessed and processed by systems such as Longleaf

For inactive data, please move it to /proj - contact us if you have no place or insufficient quota to move inactive data off of /work.

File systems such as /work in an academic research setting typically employ a file deletion policy, auto-deleting files of a certain age. **At this time, there are no time limits for files on /work.** (but still delete your old stuff)

### Modules
See [here](https://help.rc.unc.edu/modules)
