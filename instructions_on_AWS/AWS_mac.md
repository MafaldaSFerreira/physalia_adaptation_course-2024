If you are using a Mac or Linux machine, you will need to open a terminal window and then use the `ssh` command. The command *ssh* stands for *secure shell* and is a way of connecting to and interacting with remote servers. You will need to log in to the cluster using both *ssh* and a *keyfile* that has been generated for you.

Firstly, download the keyfile and open a terminal window. Then copy the keyfile into your home directory like so:
```bash
cp anyname.pem ~
```
Then, you should be able to log in with ssh. You need to provide ssh with the path to your key, which you can do with the `-i` flag. This basically points to your identity file or keyfile. For example:
```bash
ssh -i "~/anyname.pem" anyname@54.245.175.86
```
**Note that you will need to change the log in credentials shown here (i.e. the username and keyfile name) with your own**. Also be aware that the cluster IP address will change everyday. We will update you on this each day. You might be prompted to accept an RSA key - if so, just type yes and you will log in to the cluster!

### Downloading and uploading files
Occassionally, we will need to transfer files between the cluster computer and our local computer. To do this, we can use a command utility called `scp`, which stans for *secure copy* and has this syntax: 
```bash
scp [options] [[user@]host1:]source_file_or_directory [[user@]host2:]destination
```
It uses the same authentication and security as with `ssh`.

Let’s make a dummy file in our local home directory and then upload it to our home directory on the cluster using `scp`:
```bash
# make a file
touch test_file
# upload to cluster
scp -i "~/anyname.pem" test_file anyname@54.245.175.86:~/
```
Here we provide the path to the key (`-i "~/anyname.pem"`). Then, indicate where the local source file is (`test_file`), provide information to connect to the cluster and after the `:` symbol, we  specify the destination on the cluster (where we want the file to be located). Here we use `~/` to specify the home directory (`anyname@54.245.175.86:~/`).

Copying files from the cluster to our local machine is just as straightforward. You can do that like so:
```bash
# download to local
scp -i "~/mark.pem" mark@54.245.175.86:~/test_file ./
```
Where here all we did was use scp with the cluster address and path first, and the location on our computer (in our working directory) second (i.e. `./`).

Alternatively, we can use software with a graphical user interface such as Cyberduck, FileZilla, or similar, to transfer data back and forth.

Open Cyberduck, and click on 'Open connection' on the top left of the screen. Select 'SFTP (SSH File Tranfer Protocol)' from the dropdown menu and copy the IP address for the day in the Server field and your username for in the Username field. Leave the Password field empty, go down to SSH Private Key and add the Private Key Carlo sent you. Press Connect and you're connected!
