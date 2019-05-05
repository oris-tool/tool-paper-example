### A ready-to-use project on a software rejuvenation example

This repository provides a ready-to-use Maven project that you can easily import into an Eclipse workspace to experiment with the model of software rejuvenation discussed in the following paper:
"The ORIS Tool: Quantitative Evaluation of Non-Markovian Systems",
by Marco Paolieri, Marco Biagi, Laura Carnevali, Enrico Vicario,
IEEE Transaction on Software Engineering,
to appear.

Just follow these steps:

1. **Install Java 9.** For Windows and macOS, you can download a
   [package from Oracle](http://www.oracle.com/technetwork/java/javase/downloads/jdk9-downloads-3848520.html); on Debian unstable (sid) or testing (buster), or Ubuntu "bionic", you can just run `apt-get install openjdk-9-jdk`.

2. **Download Eclipse.** The [Eclipse IDE for Java Developers](http://www.eclipse.org/downloads/eclipse-packages/) package is sufficient.

3. **Clone this project.** Inside Eclipse:
   - Select `File > Import > Maven > Check out Maven Projects from
     SCM` and click `Next`.
   - If the `SCM URL` dropbox is grayed out, click on `m2e
     Marketplace` and install `m2e-egit`. You will have to restart
     Eclipse (be patient...).
   - As `SCM URL`, type:
     `https://github.com/oris-tool/tool-paper-example.git` and click
     `Next` and then `Finish`.

Your Eclipse project is ready! Just navigate to `src/main/java` and open `SoftwareRejuvenation.java` inside the package `org.oristool.examples`.
