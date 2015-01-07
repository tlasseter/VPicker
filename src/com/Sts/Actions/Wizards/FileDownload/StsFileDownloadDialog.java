package com.Sts.Actions.Wizards.FileDownload;

import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;

/*
 * ProgressBarDemo.java is a 1.4 application that requires these files:
 *   LongTask.java
 *   SwingWorker.java
 */
public class StsFileDownloadDialog extends JPanel implements ActionListener {
    public final static int ONE_SECOND = 1000;

    private String currentDirectory = null;
    private JFileChooser chooseDirectory = null;
    private JFrame frame = null;
    StsJar jar = null;
    private StsCursor cursor = null;

    int progress = 0;

    protected String ofilename = "StsStratton.sgy";
    protected OutputStream os = null;
    protected BufferedOutputStream bos = null;
    protected ByteArrayOutputStream baos = null;
    protected DataOutputStream dos = null;

    String ifilename = "S2SStratton.jar";
    int filesize = 124883600;
    static protected final int bufSize = 4096000;

    protected InputStream is = null;
    protected BufferedInputStream bis = null;
    protected DataInputStream dis = null;

    private Timer timer;
    private processFile task;

    private JLabel taskOutput = new JLabel("");

    JPanel jPanel1 = new JPanel();
    JPanel jPanel2 = new JPanel();
    JPanel mainPanel = new JPanel();
    JLabel title = new JLabel();
    JTextField segyDirectoryTxt = new JTextField();
    JTextArea jTextArea1 = new JTextArea();
    JButton directoryBrowseButton = new JButton();
    JButton acceptBtn = new JButton();
    JButton cancelBtn = new JButton();
    JButton exitBtn = new JButton();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    GridBagLayout gridBagLayout3 = new GridBagLayout();
    JProgressBar progressBar;
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout4 = new GridBagLayout();

    private String newline = "\n";

    public StsFileDownloadDialog()
    {
        super(new BorderLayout());
        task = new processFile();

        progressBar = new JProgressBar(0, task.getLengthOfTask());
        progressBar.setValue(0);
        progressBar.setStringPainted(true);
        progressBar.setSize(new Dimension(15,200));

        this.setLocation(200,200);
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setLayout(gridBagLayout4);
        jPanel1.setBackground(SystemColor.menu);
        jPanel2.setBorder(BorderFactory.createEtchedBorder());
        jPanel2.setLayout(gridBagLayout2);
        jPanel2.setBackground(SystemColor.menu);
        title.setFont(new java.awt.Font("Dialog", 1, 14));
        title.setForeground(new Color(0, 55, 152));
        title.setOpaque(false);
        title.setHorizontalAlignment(SwingConstants.CENTER);
        title.setText("Stratton Field SEGY Dataset");
        segyDirectoryTxt.setText("directory name");
        segyDirectoryTxt.setBackground(Color.WHITE);
        jTextArea1.setBackground(Color.lightGray);
        jTextArea1.setFont(new java.awt.Font("Dialog", 0, 10));
        jTextArea1.setToolTipText("");
        jTextArea1.setText("Select the directory where the SegY file will be copied. A minimum " +
                           "of 130 Mbytes is required.");
        jTextArea1.setLineWrap(true);
        jTextArea1.setWrapStyleWord(true);
        ImageIcon icon = StsIcon.createIcon("dir16x32.gif");
        directoryBrowseButton.setIcon(icon);
        directoryBrowseButton.setBackground(SystemColor.menu);
        directoryBrowseButton.setBorder(BorderFactory.createRaisedBevelBorder());

        acceptBtn.setText("Accept");
        cancelBtn.setText("Cancel");
        exitBtn.setText("Exit");
        acceptBtn.setBackground(SystemColor.menu);
        cancelBtn.setBackground(SystemColor.menu);
        exitBtn.setBackground(SystemColor.menu);

        mainPanel.setLayout(gridBagLayout1);
        mainPanel.setBackground(SystemColor.menu);
        mainPanel.setMinimumSize(new Dimension(250, 300));
        mainPanel.setPreferredSize(new Dimension(400, 200));

        mainPanel.add(jPanel1,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 0, 0), 5, 6));
        jPanel1.add(title,  new GridBagConstraints(0, 0, 2, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 4, 0, 5), 181, 11));
        jPanel1.add(jTextArea1,  new GridBagConstraints(0, 1, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 17, 0, 14), 0, 6));
        jPanel1.add(segyDirectoryTxt,  new GridBagConstraints(1, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(6, 8, 0, 5), 224, 0));
        jPanel1.add(directoryBrowseButton,  new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(7, 19, 0, 0), 10, -3));
        jPanel1.add(progressBar,   new GridBagConstraints(0, 3, 2, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(10, 16, 0, 5), 240, 0));
        jPanel1.add(taskOutput,   new GridBagConstraints(0, 4, 2, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 16, 5, 5), 240, 0));
        mainPanel.add(jPanel2,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 5, 5, 0), 250, 0));
        jPanel2.add(acceptBtn, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 0, 0));
        jPanel2.add(cancelBtn,   new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 0, 5, 5), 0, 0));
//        jPanel2.add(exitBtn,   new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
//            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 0, 5, 5), 0, 0));
        add(mainPanel, BorderLayout.WEST);

        directoryBrowseButton.addActionListener(this);
        acceptBtn.addActionListener(this);
        cancelBtn.addActionListener(this);

        //Create a timer.
        timer = new Timer(ONE_SECOND, new ActionListener()
        {
            public void actionPerformed(ActionEvent evt) {
                progressBar.setValue(task.getCurrent());
                String s = task.getMessage();
                if(s != null)
                {
                    taskOutput.setText(s);
                }
                if(task.isDone()) {
                    Toolkit.getDefaultToolkit().beep();
                    timer.stop();
                    acceptBtn.setEnabled(true);
                    setCursor(null); //turn off the wait cursor
                    progressBar.setValue(progressBar.getMinimum());
                }
            }
        });

    }

    /**
     * Called when the user presses the start button.
     */
    public void actionPerformed(ActionEvent evt) {
        Object source = evt.getSource();
        if (source == directoryBrowseButton)
        {
            System.out.println("Choose Directory...");
            if (chooseDirectory == null)
            {
                chooseDirectory = new JFileChooser(currentDirectory);
                chooseDirectory.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            }

            chooseDirectory = new JFileChooser(currentDirectory);
            chooseDirectory.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            chooseDirectory.setDialogTitle("Select or Enter Desired Directory and Press Open");
            chooseDirectory.setApproveButtonText("Open Directory");
            while(true)
            {
                chooseDirectory.showOpenDialog(null);
                if(chooseDirectory.getSelectedFile() == null)
                    break;
                File newDirectory = chooseDirectory.getSelectedFile();
                if (newDirectory.isDirectory())
                {
                    setCurrentDirectory(newDirectory.getAbsolutePath());
                    break;
                }
                else
                {
                    // File or nothing selected, strip off file name
                    String dirString = newDirectory.getPath();
                    newDirectory = new File(dirString.substring(0,dirString.lastIndexOf(File.separator)));
                    if(newDirectory.isDirectory())
                    {
                        setCurrentDirectory(newDirectory.getAbsolutePath());
                        break;
                    }
                    if(!StsYesNoDialog.questionValue(this,"Must select the directory that\n contains the SegY Files.\n\n Continue?"))
                        break;
                }
            }
        }
        else if(source == acceptBtn)
        {
            acceptBtn.setEnabled(false);
            setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
            task.go();
            timer.start();
        }
        else if(source == cancelBtn)
        {
            acceptBtn.setEnabled(true);
            setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
            task.stop();
            progressBar.setValue(0);
            taskOutput.setText("");
        }
        else if(source == exitBtn)
        {
            frame.setVisible(false);
        }
    }

    private void setCurrentDirectory(String directory)
    {
        currentDirectory = directory;
        segyDirectoryTxt.setText(currentDirectory);
    }

    /**
     * Create the GUI and show it.  For thread safety,
     * this method should be invoked from the
     * event-dispatching thread.
     */
    private static void createAndShowGUI() {
        //Make sure we have nice window decorations.
//        JFrame.setDefaultLookAndFeelDecorated(true);

        //Create and set up the window.
        JFrame frame = new JFrame("Segy File Download");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        //Create and set up the content pane.
        JComponent newContentPane = new StsFileDownloadDialog();
        newContentPane.setOpaque(true); //content panes must be opaque
        frame.setContentPane(newContentPane);

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }

    public static void main(String[] args)
    {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { createAndShowGUI(); }});
    }

private class processFile {
    private int lengthOfTask;
    private int current = 0;
    private boolean done = false;
    private boolean canceled = false;
    private String statMessage;

    public processFile() {
        lengthOfTask = filesize;
    }

    /**
     * Called from ProgressBarDemo to start the task.
     */
    public void go() {
        final SwingWorker worker = new SwingWorker() {
            public Object construct() {
                current = 0;
                done = false;
                canceled = false;
                statMessage = null;
                return new extractFile();
            }
        };
        worker.start();
    }

    /**
     * Called from ProgressBarDemo to find out how much work needs
     * to be done.
     */
    public int getLengthOfTask() {
        return lengthOfTask;
    }

    /**
     * Called from ProgressBarDemo to find out how much has been done.
     */
    public int getCurrent() {
        return current;
    }

    public void stop() {
        current = 0;
        canceled = true;
        done = true;
        statMessage = null;
    }

    /**
     * Called from ProgressBarDemo to find out if the task has completed.
     */
    public boolean isDone() {
        return done;
    }

    /**
     * Returns the most recent status message, or null
     * if there is no current status message.
     */
    public String getMessage() {
        return statMessage;
    }

    /**
     * The actual long running task.  This runs in a SwingWorker thread.
     */
    class extractFile
    {
        extractFile()
        {
            jar = (StsJar) StsWebStartJar.constructor(ifilename);
            if (jar == null)
            {
                statMessage = "Failed to find SegY jar file: S2SStratton.jar in cache, try download again.";
                return;
            }
            String[] filenames = jar.getFilenames();
            if (filenames == null || filenames.length < 1)
            {
                statMessage = "Failed to find Stratton SegY in jar file. File might be corrupt, try download again.";
                return;
            }
            else
            {
                StsAbstractFile iFile = jar.getFile(0);
                System.out.println("Input file=" + iFile.getFilename());
                StsFile oFile = StsFile.constructor(currentDirectory +
                    File.separator + ofilename);
                System.out.println("Output file:" + oFile.getFilename());

                //
                // Read data from jar and output to user directory
                //
                try
                {
                    is = iFile.getInputStream();
                }
                catch(Exception e)
                {
                    new StsMessage(StsModel.getCurrentModel().win3d,  StsMessage.WARNING, "Failed to find file " + iFile.getFilename());
                    return;
                }
                bis = new BufferedInputStream(is, bufSize);
                dis = new DataInputStream(bis);

                try
                {
                    os = oFile.getOutputStream();
                }
                catch(Exception e)
                {
                    new StsMessage(StsModel.getCurrentModel().win3d,  StsMessage.WARNING, "Failed to find file " + oFile.getFilename());
                    return;
                }
                bos = new BufferedOutputStream(os, bufSize);
                dos = new DataOutputStream(bos);

                byte[] buffer = new byte[bufSize];
                int read = 0;
                try
                {
                    while(!canceled && !done)
                    {
                        read = dis.read(buffer);

                        if(read <= 0)
                        {
                            done = true;
                            current = lengthOfTask;
                            break;
                        }
                        dos.write(buffer, 0, read);
                        if(!canceled)
                        {
                            current += read;
                            statMessage = "Completed " + current + " out of " + lengthOfTask + " bytes.";
                        }
                    }
                    dos.flush();
                }
                catch (Exception e) {
                    statMessage = "Failed to copy the Stratton SegY file. File might be corrupt, try download again.";
                    return;
                }
            }
        }
    }
}


/**
 * This is the 3rd version of SwingWorker (also known as
 * SwingWorker 3), an abstract class that you subclass to
 * perform GUI-related work in a dedicated thread.  For
 * instructions on and examples of using this class, see:
 *
 * http://java.sun.com/docs/books/tutorial/uiswing/misc/threads.html
 *
 * Note that the API changed slightly in the 3rd version:
 * You must now invoke start() on the SwingWorker after
 * creating it.
 */
private abstract class SwingWorker {
    private Object value;  // see getValue(), setValue()

    /**
     * Class to maintain reference to current worker thread
     * under separate synchronization control.
     */
    private class ThreadVar {
        private Thread thread;
        ThreadVar(Thread t) { thread = t; }
        synchronized Thread get() { return thread; }
        synchronized void clear() { thread = null; }
    }

    private ThreadVar threadVar;

    /**
     * Get the value produced by the worker thread, or null if it
     * hasn't been constructed yet.
     */
    protected synchronized Object getValue() {
        return value;
    }

    /**
     * Set the value produced by worker thread
     */
    private synchronized void setValue(Object x) {
        value = x;
    }

    /**
     * Compute the value to be returned by the <code>get</code> method.
     */
    public abstract Object construct();

    /**
     * Called on the event dispatching thread (not on the worker thread)
     * after the <code>construct</code> method has returned.
     */
    public void finished() {
    }

    /**
     * A new method that interrupts the worker thread.  Call this method
     * to force the worker to stop what it's doing.
     */
    public void interrupt() {
        Thread t = threadVar.get();
        if (t != null) {
            t.interrupt();
        }
        threadVar.clear();
    }

    /**
     * Return the value created by the <code>construct</code> method.
     * Returns null if either the constructing thread or the current
     * thread was interrupted before a value was produced.
     *
     * @return the value created by the <code>construct</code> method
     */
    public Object get() {
        while (true) {
            Thread t = threadVar.get();
            if (t == null) {
                return getValue();
            }
            try {
                t.join();
            }
            catch (InterruptedException e) {
                Thread.currentThread().interrupt(); // propagate
                return null;
            }
        }
    }


    /**
     * Start a thread that will call the <code>construct</code> method
     * and then exit.
     */
    public SwingWorker()
    {
        final Runnable doFinished = new Runnable()
        {
           public void run() { finished(); }
        };

        Runnable doConstruct = new Runnable()
        {
            public void run()
            {
                try
                {
                    setValue(construct());
                }
                finally
                {
                    threadVar.clear();
                }

                SwingUtilities.invokeLater(doFinished);
            }
        };

        Thread t = new Thread(doConstruct);
        threadVar = new ThreadVar(t);
    }

    /**
     * Start the worker thread.
     */
    public void start() {
        Thread t = threadVar.get();
        if (t != null) {
            t.start();
        }
    }
}
}
