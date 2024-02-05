# Explaining Channels

- Nextflow is based on the dataflow programming model in which **processes communicate through channels**.
- The channel has two major properties:
    - Sending a message is an asynchronous operation
    - Recieving a message is a synchronous operation
    - Basically, a process can send a message and continue all it has to do. But the other process need to wait until the message has arrived to continue. 
- The declaration of a channel can be before the workflow scope or within it. As long as it is upstream of the process that requires the specific channel (no hoisting).


