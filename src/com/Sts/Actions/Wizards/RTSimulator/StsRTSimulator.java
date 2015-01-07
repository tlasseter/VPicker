package com.Sts.Actions.Wizards.RTSimulator;

import com.Sts.Actions.Wizards.SensorLoad.StsSensorFile;
import com.Sts.DBTypes.StsMonitor;
import com.Sts.IO.StsAbstractFile;
import com.Sts.IO.StsAsciiFile;
import com.Sts.IO.StsFile;
import com.Sts.Interfaces.StsProgressRunnable;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Icons.StsIcon;
import com.Sts.UI.Progress.StsProgressPanel;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsYesNoDialog;
import com.Sts.Utilities.DateTime.CalendarParser;
import com.Sts.Utilities.StsToolkit;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.Date;

/*
 * ProgressBarDemo.java is a 1.4 application that requires these files:
 *   LongTask.java
 *   SwingWorker.java
 */

public class StsRTSimulator {

    /**
     * Create the GUI and show it.  For thread safety,
     * this method should be invoked from the
     * event-dispatching thread.
     */
    private static void createAndShowGUI() {
        JTabbedPane tabbedPanels = new JTabbedPane();

        JDialog dialog = new JDialog(new JFrame("Real-time Simulator"));
        StsRTSensorSimulatorPanel sensorPanel = new StsRTSensorSimulatorPanel();
        StsRTWellSimulatorPanel wellPanel = new StsRTWellSimulatorPanel();

        tabbedPanels.add("Sensors", sensorPanel);
        tabbedPanels.add("Wells", wellPanel);

        dialog.getContentPane().add(tabbedPanels);
        dialog.pack();
        dialog.setVisible(true);
    }

    public static void main(String[] args)
    {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        StsToolkit.runLaterOnEventThread(new Runnable()
        {
            public void run()
            { createAndShowGUI(); }
        });
    }

}