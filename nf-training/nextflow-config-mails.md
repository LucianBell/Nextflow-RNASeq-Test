# Mail configuration

- We can send a notification email when the workflow execution completes using the `-N <email address>` **command-line option**.
- However, **this requires the configuration of a SMTP server in the nextflow config file.**
- You have an example on how to do that on the `nextflow.config` file

## Workflow function
- The built-in function `sendMail` allows you to send a mail message from a workflow script.
- Example:

        sendMail(
            to: 'you@gmail.com',
            subject: 'Catch up',
            body: 'Hi, how are you!',
            attach: '/some/path/attachment/file.txt'
        )

- Check more about this function [here](https://www.nextflow.io/docs/latest/mail.html#mail-configuration)