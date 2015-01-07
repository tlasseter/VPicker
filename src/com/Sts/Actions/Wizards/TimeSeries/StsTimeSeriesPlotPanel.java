package com.Sts.Actions.Wizards.TimeSeries;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import info.monitorenter.gui.chart.*;
import info.monitorenter.gui.chart.controls.*;
import info.monitorenter.gui.chart.labelformatters.*;
import info.monitorenter.gui.chart.traces.*;
import info.monitorenter.gui.chart.views.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.text.*;
import java.util.*;
/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */
public class StsTimeSeriesPlotPanel extends JPanel implements MouseListener, MouseMotionListener
{
    private StsTimeSeriesWizard wizard;
    private StsTimeSeriesPlot wizardStep;

    private Chart2D chart = null;
    private ChartPanel chartPanel = null;
    private ITrace2D[] traces = null;

    private Graphics gc = null;

    private float[][] attributes = null;
    private String[] names = null;
    private float[] time = null;
    private String[] timeString = null;
    boolean debug = false;
    boolean movieInitialized = false;

    int movieInc = 600000;
    long movieStart = 0;
    StsMovie movie = null;
    int nFrames = 0;
    int currentMovieFrame = 0;
    int startFrame = 0, lastFrame = 0;;
    boolean enableDraw = true;
    int frameInc = 1;
    long startMouseY = -1;
    long endMouseY = -1;
    boolean dragging = false;

    JPanel jPanel1 = new JPanel();
    JLabel timeLabel = new JLabel("Time: 0");

    SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd kk:mm:ss.S");

    GridBagLayout gridBagLayout2 = new GridBagLayout();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    BorderLayout borderLayout1 = new BorderLayout();

    int currentItem = 0, startItem = 0, endItem = 0;

    public StsTimeSeriesPlotPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsTimeSeriesWizard)wizard;
        this.wizardStep = (StsTimeSeriesPlot)wizardStep;

        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        getAttributeData();

        traces = new Trace2DSimple[wizard.getAttributeIndices().length];
        chart = createChart();

        new LayoutFactory.BasicPropertyAdaptSupport(this.getRootPane(), chart);
        LayoutFactory factory = LayoutFactory.getInstance();
        factory.setShowAxisXRangePolicyMenu(false);
        factory.setShowSaveImageMenu(false);
        factory.setShowAxisYRangePolicyMenu(false);
        factory.setShowAxisXRangePolicyMenu(false);

        factory.setShowTraceZindexMenu(false);
        factory.setShowTraceSelectionMenu(true);

        chartPanel = new ChartPanel(chart);

        jPanel1.add(chartPanel,   new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));

        wizard.rebuild();
        chart.setForeground(Color.gray);
        chart.addMouseListener(this);
        chart.addMouseMotionListener(this);
    }

    public void destroyChart()
    {
        time = null;
        timeString = null;
        chart = null;
        traces = null;
        attributes = null;
        currentItem = 0;
        movieInitialized = false;
        repaint();
    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout1);

        jPanel1.setLayout(gridBagLayout2);
        jPanel1.setBorder(BorderFactory.createEtchedBorder());

        timeLabel.setFont(new Font("Dialog", Font.BOLD, 14));
        this.add(timeLabel,   new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));
        this.add(jPanel1,   new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));
    }

    private Chart2D createChart()
    {
        Chart2D chart1 = new Chart2D();
        for(int i=0; i<names.length; i++)
        {
            Color color = StsColor.getColorByName(StsColor.colorNames[i+1]).getColor();
            traces[i] = createDataset(names[i], attributes[i]);
            traces[i].setColor(color);
            chart1.addTrace(traces[i]);
            IAxis axis = chart1.getAxisY();
            axis.setPaintGrid(true);
            axis.setFormatter(new LabelFormatterDate(new SimpleDateFormat("hh:mm:ss.S")));
            axis = chart1.getAxisX();
            axis.setPaintGrid(true);
        }
        chart1.getCurrentTrace().setSelected(true);
        chart1.setBackground(Color.BLACK);
        chart1.setForeground(Color.WHITE);
        chart1.setPaintLabels(true);
        chart1.setGridColor(Color.DARK_GRAY);
        chart1.setPaintCrosshair(true);
        return chart1;
    }

    private Trace2DSimple createDataset(String name, float[] attribute)
    {

        Date date = new Date(0L);
        Trace2DSimple trace = new Trace2DSimple(name);

        try
        {
            if(debug) System.out.println("Number of attributes:" + attribute.length);
            long startTime = wizard.getAsciiStartTime();
            for (int i = 0; i < attribute.length; i++)
            {
                if(wizard.getTimeType() == wizard.ARTIFICIAL_TIME)
                {
                    date = new Date( (long) (startTime + (long) (time[i] * 60 * 1000)));
                }
                else
                {
                    if(timeString[i] == null)
                        break;
                    timeString[i] = wizard.cleanTimeString(timeString[i]);
                    if(debug) System.out.println("TimeString=" + timeString[i]);
                    if(timeString[i] == null)
                    {
                        StsMessageFiles.infoMessage("Invalid time value (" + time + ") in time series file, skipping value");
                        continue;
                    }
                    date = format.parse(timeString[i]);
                }
                if(debug) System.out.println("Date[" + i + "]=" + date.getTime() + " Attribute[" + i + "]=" + attribute[i]);
                trace.addPoint(attribute[i], date.getTime());
            }
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage("Failed to create date, setting to 0L");
            return null;
        }
        return trace;
    }

    public boolean initializeMovie()
    {
        movie = wizard.getMovie();
        movieInc = wizard.getMovieIncrement() * (int)movie.getIncrement();

        nFrames = movie.getNumberFrames();
        if(wizard.isMovieStartClockTime())
            movieStart = wizard.getMovieStart();
        else
            movieStart = (long)traces[0].getMinY() + wizard.getMovieStart();
        long movieEnd = movieStart + (movieInc * (nFrames-1));

        movie.initializeEnvironment();
        currentMovieFrame = 0;
        startItem = 0;
        endItem = traces[0].getSize();

        long start2d = (long)traces[0].getMinY();
        long end2d = (long)traces[0].getMaxY();
        while(movieStart < start2d)
        {
            movieStart = movieStart + movieInc;
            startFrame++;
        }
        lastFrame = movie.getNumberFrames() - 1;
        while(movieEnd > end2d)
        {
            movieEnd = movieEnd - movieInc;
            lastFrame--;
        }
        if(debug) System.out.println("StartFrame=" + startFrame + " LastFrame=" + lastFrame);
        movieInitialized = true;
        startFrame();
        return true;
    }

    public void setTimeLabel(long time)
    {
        Date date = new Date(time);
        String timeLbl = new SimpleDateFormat("hh:mm:ss.S").format(date);
        timeLabel.setText(timeLbl);
    }

    public void nextFrame()
    {
        chart.clearXColorRange();
        long time = computeCurrentItem();
        setTimeLabel(time);
        if((currentMovieFrame >= startFrame) && (currentMovieFrame <= lastFrame))
            chart.setCrosshairXValue(time);
        if(movie.getElevationAnimation())
            movie.setElevationAndAzimuth(movie.INCREMENT, movie.ELEVATION);
        if(movie.getAzimuthAnimation())
            movie.setElevationAndAzimuth(movie.INCREMENT, movie.AZIMUTH);

        movie.setCurrentFrame(currentMovieFrame);

        if(currentMovieFrame == movie.getNumberFrames()-1)
            currentMovieFrame = 0;
        else
            currentMovieFrame++;

    }

    private long computeCurrentItem()
    {
        if(wizard.isMovieStartClockTime())
            movieStart = wizard.getMovieStart();
        else
            movieStart = (long)traces[0].getMinY() + wizard.getMovieStart();

        long millis = movieStart + (currentMovieFrame * movieInc);
//        currentItem = series[0].index(new Millisecond(new Date(millis)));   NEED TO COMPUTE CURRENT VALUE
        currentItem = 0;
        if(currentItem != -1)
            currentItem = Math.abs(currentItem);

        if(debug) System.out.println("CurrentItem= " + currentItem);
        return millis;
    }

    public void previousFrame()
    {
        chart.clearXColorRange();
        long time = computeCurrentItem();
        setTimeLabel(time);

        if((currentMovieFrame >= startFrame) && (currentMovieFrame <= lastFrame))
            chart.setCrosshairXValue(time);

        if(movie.getElevationAnimation())
            movie.setElevationAndAzimuth(movie.INCREMENT, movie.ELEVATION);
        if(movie.getAzimuthAnimation())
            movie.setElevationAndAzimuth(movie.INCREMENT, movie.AZIMUTH);

       movie.setCurrentFrame(currentMovieFrame);

        if(currentMovieFrame == 0)
            currentMovieFrame = movie.getNumberFrames() - 1;
        else
            currentMovieFrame--;
    }

    public void stop()
    {
        chart.clearXColorRange();
        currentMovieFrame = 0;
        long time = computeCurrentItem();
        setTimeLabel(time);

        if((currentMovieFrame >= startFrame) && (currentMovieFrame <= lastFrame))
            chart.setCrosshairXValue(time);
        movie.setCurrentFrame(currentMovieFrame);
    }

    public void startFrame()
    {
        chart.clearXColorRange();
        currentMovieFrame = 0;
        long time = computeCurrentItem();
        setTimeLabel(time);

        if((currentMovieFrame >= startFrame) && (currentMovieFrame <= lastFrame))
            chart.setCrosshairXValue(time);
        movie.setCurrentFrame(currentMovieFrame);
    }

    public void endFrame()
    {
        chart.clearXColorRange();
        currentMovieFrame = movie.getNumberFrames()-1;
        long time = computeCurrentItem();
        setTimeLabel(time);

        if((currentMovieFrame >= startFrame) && (currentMovieFrame <= lastFrame))
            chart.setCrosshairXValue(time);

        movie.setCurrentFrame(currentMovieFrame);
    }

    private void getAttributeData()
    {
        BufferedReader bufRdr = null;
        String line = null;
        StringTokenizer stok = null;
        boolean firstRow = true;
        File asciiFile = null;
        int attIdx = 0, nAtts=0;
        int[] attIndices = null;

        attIndices = wizard.getAttributeIndices();
        nAtts = attIndices.length;
        attributes = new float[nAtts][wizard.getNumberValidRows()];
        if(wizard.getTimeType() == wizard.ARTIFICIAL_TIME)
            time = new float[wizard.getNumberValidRows()];
        else
            timeString = new String[wizard.getNumberValidRows()];
        names = new String[nAtts];

        try
        {
            asciiFile = wizard.getAsciiFile();
            bufRdr = new BufferedReader(new FileReader(asciiFile));
            int nRows = 0;
            while((line = bufRdr.readLine()) != null)
            {
                stok = new StringTokenizer(line,", ;");
                int nTokens = stok.countTokens();
                if(nTokens != wizard.getNumberTokens())
                    continue;

                // Create column labels if exist in file
                if(firstRow)
                {
                    firstRow = false;
                    attIdx = 0;
                    for(int i=0; i<wizard.getNumberTokens(); i++)
                    {
                        String token = stok.nextToken();
                        if(attIndices[attIdx] == i)
                        {
                            names[attIdx] = token;
                            attIdx++;
                        }
                        if(attIdx >= attIndices.length)
                            break;
                    }
                    continue;
                }
                attIdx = 0;
                for(int i=0; i<wizard.getNumberTokens(); i++)
                {
                    String token = stok.nextToken();
                    if(attIdx < attIndices.length)
                    {
                        if (attIndices[attIdx] == i)
                        {
                            attributes[attIdx][nRows] = Float.parseFloat(token);
                            attIdx++;
                        }
                    }
                    if(wizard.getTimeIndex() == i)
                    {
                        if(wizard.getTimeType() == wizard.ARTIFICIAL_TIME)
                            time[nRows] = Float.parseFloat(token);
                        else
                        {
                            timeString[nRows] = token;
                        }
                    }
                }
                nRows++;
                if(nRows >= wizard.getNumberValidRows() - 1)
                    break;
            }
            bufRdr.close();
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage("Failed to read file: " + asciiFile.getName());
            return;
        }
    }

    public void mouseExited(MouseEvent e) { }
    public void mouseEntered(MouseEvent e) { }
    public void mouseReleased(MouseEvent e)
    {
        startMouseY = -1;
        if (chart != null)
        {
            endMouseY = (long)chart.getAxisY().translatePxToValue(e.getY());
            if(debug)
                System.out.println("Selected time=" + endMouseY);
            if(debug)
                System.out.println("Date=" + new Date(endMouseY).toString());

            currentMovieFrame = (int)((long)(endMouseY - movieStart)/movieInc);

            if(currentMovieFrame < 1)
                currentMovieFrame = 0;
            if(currentMovieFrame > nFrames)
                currentMovieFrame = movie.getNumberFrames() - 1;
            if(debug)
                System.out.println("User button in graph event --> Current frame= " + currentMovieFrame);

            if(!dragging)
                movie.setCurrentFrame(currentMovieFrame);
            else
                dragging = false;

            chart.setCrosshairXValue(endMouseY);
            computeCurrentItem();
        }
    }
    public void mouseClicked(MouseEvent e)
    {
    }
    public void mousePressed(MouseEvent e)
    {
        chart.clearXColorRange();
        startMouseY = (long)chart.getAxisY().translatePxToValue(e.getY());
        chart.setCrosshairXValue(startMouseY);
    }
    public void mouseDragged(MouseEvent e)
    {
        if((chart != null) && (startMouseY != -1))
        {
            dragging = true;
            endMouseY = (long)chart.getAxisY().translatePxToValue(e.getY());
            currentMovieFrame = (int)((long)(endMouseY - movieStart)/movieInc);

            if(currentMovieFrame < 1)
                currentMovieFrame = 0;
            if(currentMovieFrame > nFrames)
                currentMovieFrame = movie.getNumberFrames() - 1;

            if(endMouseY < startMouseY)
            {
                chart.setXColorRange(startMouseY, endMouseY);
                movie.setFrameRange(currentMovieFrame, (int)((long)(startMouseY - movieStart)/movieInc));
            }
            if(startMouseY < endMouseY)
            {
                chart.setXColorRange(endMouseY, startMouseY);
                movie.setFrameRange((int)((long)(startMouseY - movieStart)/movieInc), currentMovieFrame);
            }
        }
    }
    public void mouseMoved(MouseEvent e)
    {

    }
}
