nextflow.enable.dsl=2

process sayHello {
    output:
    stdout

    """
    echo 'Hello from Nextflow!'
    """
}

workflow {
    sayHello()
}

