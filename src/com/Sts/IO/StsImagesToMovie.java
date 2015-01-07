package com.Sts.IO;
/*
 * @(#)JpegImagesToMovie.java	1.3 01/03/13
 *
 * Copyright (c) 1999-2001 Sun Microsystems, Inc. All Rights Reserved.
 *
 * Sun grants you ("Licensee") a non-exclusive, royalty free, license to use,
 * modify and redistribute this software in source and binary code form,
 * provided that i) this copyright notice and license appear on all copies of
 * the software; and ii) Licensee does not utilize the software in a manner
 * which is disparaging to Sun.
 *
 * This software is provided "AS IS," without a warranty of any kind. ALL
 * EXPRESS OR IMPLIED CONDITIONS, REPRESENTATIONS AND WARRANTIES, INCLUDING ANY
 * IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE OR
 * NON-INFRINGEMENT, ARE HEREBY EXCLUDED. SUN AND ITS LICENSORS SHALL NOT BE
 * LIABLE FOR ANY DAMAGES SUFFERED BY LICENSEE AS A RESULT OF USING, MODIFYING
 * OR DISTRIBUTING THE SOFTWARE OR ITS DERIVATIVES. IN NO EVENT WILL SUN OR ITS
 * LICENSORS BE LIABLE FOR ANY LOST REVENUE, PROFIT OR DATA, OR FOR DIRECT,
 * INDIRECT, SPECIAL, CONSEQUENTIAL, INCIDENTAL OR PUNITIVE DAMAGES, HOWEVER
 * CAUSED AND REGARDLESS OF THE THEORY OF LIABILITY, ARISING OUT OF THE USE OF
 * OR INABILITY TO USE SOFTWARE, EVEN IF SUN HAS BEEN ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGES.
 *
 * This software is not designed or intended for use in on-line control of
 * aircraft, air traffic, aircraft navigation or aircraft communications; or in
 * the design, construction, operation or maintenance of any nuclear
 * facility. Licensee represents and warrants that it will not use or
 * redistribute the Software for such purposes.
 */

import com.Sts.UI.Progress.*;

import javax.media.*;
import javax.media.control.*;
import javax.media.datasink.*;
import javax.media.format.*;
import javax.media.protocol.*;
import java.awt.*;
import java.io.*;
import java.util.*;

/**
 * This program takes a list of JPEG image files and convert them into
 * a QuickTime movie.
 */
public class StsImagesToMovie implements ControllerListener, DataSinkListener
{
    StsProgressPanel panel = null;
    static boolean debug = false;
    
    public boolean doIt(int width, int height, int frameRate, Vector inFiles, MediaLocator outML, StsProgressPanel panel)
    {
        this.panel = panel;

        ImageDataSource ids = new ImageDataSource(width, height, frameRate, inFiles, panel);
        PullBufferStream[] streams = ids.getStreams();
        for(int i=0; i<streams.length; i++)
        {
        	if(debug)
        		System.out.println("Input stream[" + i + "]=" + streams[i].toString());
        }
        Processor p;
        try
        {
            if(debug) System.out.println("Creating processor");
            panel.appendLine("Create processor for the image datasource ...");
            p = Manager.createProcessor(ids);
        }
        catch (Exception e)
        {
            panel.appendLine("ERROR: Cannot create a processor from the data source.");
            return false;
        }
        if(debug) System.out.println("Adding listener to processor");
        p.addControllerListener(this);

        // Put the Processor into configured state so we can set
        // some processing options on the processor.
        if(debug) System.out.println("Configuring processor");
        p.configure();
        if (!waitForState(p, p.Configured))
        {
            panel.appendLine("ERROR: Failed to configure the processor.");
            return false;
        }

        // Set the output content descriptor to QuickTime.
        if(debug) System.out.println("Setting output to QuickTime");
        p.setContentDescriptor(new ContentDescriptor(FileTypeDescriptor.QUICKTIME));
        //p.setContentDescriptor(new ContentDescriptor(FileTypeDescriptor.MPEG));
        // Query for the processor for supported formats.
        // Then set it on the processor.
        if(debug) System.out.println("Getting track controls and formats");
        TrackControl tcs[] = p.getTrackControls();
        Format f[] = tcs[0].getSupportedFormats();
        if (f == null || f.length <= 0)
        {
            if(debug) System.out.println("Unsupportted format...");        	
            panel.appendLine("ERROR: The mux does not support the input format: " + tcs[0].getFormat() + "\n");
            return false;
        }
        tcs[0].setFormat(f[0]);

        panel.appendLine("Setting the track format to: " + f[0] + "\n");

        // We are done with programming the processor.  Let's just
        // realize it.
        p.realize();
        if(!waitForState(p, p.Realized))
        {
            panel.appendLine("ERROR: Failed to realize the processor.\n");
            return false;
        }

        // Now, we'll need to create a DataSink.
        if(debug) System.out.println("Creating datasink");
        DataSink dsink = createDataSink(p, outML);
        if(dsink == null)
        {
            if(debug) System.out.println("Failed to create datasink");
            panel.appendLine("ERROR: Failed to create a DataSink for the given output MediaLocator: " + outML);
            return false;
        }
        if(debug) System.out.println("Adding listener to datasink");
        dsink.addDataSinkListener(this);
        fileDone = false;

        panel.appendLine("Starting file processing...");

        // OK, we can now start the actual transcoding.
        try
        {
            if(debug) System.out.println("Processing input files..");
            p.start();
            dsink.start();
        }
        catch (IOException e)
        {
            if(debug) System.out.println("Error processing input files");
            panel.appendLine("ERROR: IO error during processing");
            return false;
        }

        // Wait for EndOfStream event.
        waitForFileDone();

        // Cleanup.
        try
        {
            if(debug) System.out.println("Finished processing..");        	
            dsink.close();
        }
        catch (Exception e) {}
        p.removeControllerListener(this);

        panel.appendLine("Done processing files.");
        return true;
    }


    /**
     * Create the DataSink.
     */
    DataSink createDataSink(Processor p, MediaLocator outML)
    {
        DataSource ds;
        if(debug) System.out.println("createDataSink: getDataOutput()");
        if ((ds = p.getDataOutput()) == null)
        {
            System.err.println("Something is really wrong: the processor does not have an output DataSource");
            return null;
        }

        DataSink dsink;
        try
        {
        	if(debug) System.out.println("createDataSink: Manager.createDataSink(ds, outML):\n ds=" + ds.toString() +
        			"\n outML=" + outML.toString() +
        			"\n ds.getContentType()=" + ds.getContentType());
        	
            dsink = Manager.createDataSink(ds, outML);
            
            if(debug) System.out.println("createDataSink: dsink.open():\n dsink=" + dsink.toString() + 
            		"\n dsink.getOutputLocator()" + dsink.getOutputLocator() +
            		"\n dsink.getContentType()" + dsink.getContentType() +
            		"\n dsink.getControls()" + dsink.getControls());
            dsink.open();
            if(debug) System.out.println("createDataSink: completed dsink.open():\n"); 
        }
        catch (Exception e)
        {
        	System.err.println("Error creating or opening datasink.." + e + "\nnTry local file selection.");
            panel.appendLine("Exception thrown: Cannot create the DataSink: " + e.getMessage());
            panel.setDescriptionAndLevel("Exception thrown.", StsProgressBar.ERROR);
            return null;
        }
        return dsink;
    }

    Object waitSync = new Object();
    boolean stateTransitionOK = true;

    /**
     * Block until the processor has transitioned to the given state.
     * Return false if the transition failed.
     */
    boolean waitForState(Processor p, int state)
    {
        synchronized (waitSync)
        {
            try
            {
                while (p.getState() < state && stateTransitionOK)
                    waitSync.wait();
            }
            catch (Exception e) {}
        }
        return stateTransitionOK;
    }


    /**
     * Controller Listener.
     */
    public void controllerUpdate(ControllerEvent evt)
    {
        if (evt instanceof ConfigureCompleteEvent ||
            evt instanceof RealizeCompleteEvent ||
            evt instanceof PrefetchCompleteEvent)
        {
            synchronized (waitSync)
            {
                stateTransitionOK = true;
                waitSync.notifyAll();
            }
        }
        else if (evt instanceof ResourceUnavailableEvent)
        {
            synchronized (waitSync)
            {
                stateTransitionOK = false;
                waitSync.notifyAll();
            }
        }
        else if (evt instanceof EndOfMediaEvent)
        {
            evt.getSourceController().stop();
            evt.getSourceController().close();
        }
    }

    Object waitFileSync = new Object();
    boolean fileDone = false;
    boolean fileSuccess = true;

    /**
     * Block until file writing is done.
     */
    boolean waitForFileDone()
    {
        synchronized (waitFileSync)
        {
            try
            {
                while (!fileDone)
                    waitFileSync.wait();
            }
            catch (Exception e) {}
        }
        return fileSuccess;
    }


    /**
     * Event handler for the file writer.
     */
    public void dataSinkUpdate(DataSinkEvent evt)
    {
        if (evt instanceof EndOfStreamEvent)
        {
            synchronized (waitFileSync)
            {
                fileDone = true;
                waitFileSync.notifyAll();
            }
        }
        else if (evt instanceof DataSinkErrorEvent)
        {
            synchronized (waitFileSync)
            {
                fileDone = true;
                fileSuccess = false;
                waitFileSync.notifyAll();
            }
        }
    }

    public static boolean imagesToMovie(int wide, int high, int fr, String output, String filenames[], StsProgressPanel panel)
    {
        int i = 0;
        int width = -1, height = -1, frameRate = 1;
        Vector inputFiles = new Vector();
        String outputURL = null;

        width = wide;
        height = high;
        frameRate = fr;
        outputURL = output;
        /*
        System.loadLibrary("jmutil");
        System.loadLibrary("jmacm");
        System.loadLibrary("jmam");
        System.loadLibrary("jmcvid");
        System.loadLibrary("jmg723");
        System.loadLibrary("jmgdi");
        System.loadLibrary("jmgsm");
        System.loadLibrary("jmh261");
        System.loadLibrary("jmh263enc");
        System.loadLibrary("jmvh263");
        System.loadLibrary("jmjpeg");
        System.loadLibrary("jmmci");
        System.loadLibrary("jmddraw");
        System.loadLibrary("jmmpa");
        System.loadLibrary("jmmpegv");
        System.loadLibrary("jmvcm");
        System.loadLibrary("jmvfw");
        System.loadLibrary("jmdaudc");
        System.loadLibrary("jmdaud");
        System.loadLibrary("jmfjawt");
        */
        if(debug) System.out.println("Verifying output and input files.\n");        
        for(i=0; i<filenames.length; i++)
            inputFiles.addElement(filenames[i]);

        if (outputURL == null || inputFiles.size() == 0)
        {
            System.err.println("Invalid output filename or number of input files.");
            return false;
        }
        
        // Check for output file extension.
        if (!outputURL.endsWith(".mov") && !outputURL.endsWith(".MOV"))
            outputURL = outputURL + ".mov";
        //if (!outputURL.endsWith(".avi") && !outputURL.endsWith(".avi"))
        //    outputURL = outputURL + ".avi";
        if(debug) System.out.println("Building movie with name: " + outputURL);
        if (width < 0 || height < 0)
        {
            System.err.println("Please specify the correct image size.");
            width = 512;
            height = 512;
        }
        // Check the frame rate.
        if (frameRate < 1)
            frameRate = 1;

        // Generate the output media locators.
        MediaLocator oml;
        if(debug) System.out.println("Creating mediaLocator");
        if ((oml = createMediaLocator(outputURL)) == null)
        {
            System.err.println("Cannot build media locator from: " + outputURL);
            return false;
        }

        StsImagesToMovie imageToMovie = new StsImagesToMovie();
        return imageToMovie.doIt(width, height, frameRate, inputFiles, oml, panel);
    }

    /**
     * Create a media locator from the given string.
     */
    static MediaLocator createMediaLocator(String url) {

	MediaLocator ml;

	if (url.indexOf("file:") == 0 && (ml = new MediaLocator(url)) != null)
	    return ml;
	
	if (url.indexOf(":") == 1 && (ml = new MediaLocator(url)) != null)
	{
	    if ((ml = new MediaLocator("file:" + url)) != null)
		return ml;
	} else {
	    String file = "file:" + System.getProperty("user.dirNo") + File.separator + url;
	    if ((ml = new MediaLocator(file)) != null)
		return ml;
	}

	return null;
    }


    ///////////////////////////////////////////////
    //
    // Inner classes.
    ///////////////////////////////////////////////


    /**
     * A DataSource to read from a list of JPEG image files and
     * turn that into a stream of JMF buffers.
     * The DataSource is not seekable or positionable.
     */
    class ImageDataSource extends PullBufferDataSource
    {
        ImageSourceStream streams[];

        ImageDataSource(int width, int height, int frameRate, Vector images, StsProgressPanel panel)
        {
            if(debug) System.out.println("Constructing ImageSourceStream");
            streams = new ImageSourceStream[1];
            streams[0] = new ImageSourceStream(width, height, frameRate, images, panel);
        }

        public void setLocator(MediaLocator source) {}

        public MediaLocator getLocator()
        {
            return null;
        }

        /**
         * Content type is of RAW since we are sending buffers of video
         * frames without a container format.
         */
        public String getContentType()
        {
            return ContentDescriptor.RAW;
        }
        public void connect() {}
        public void disconnect() {}
        public void start() {}
        public void stop() {}

        /**
         * Return the ImageSourceStreams.
         */
        public PullBufferStream[] getStreams()
        {
            return streams;
        }

        /**
         * We could have derived the duration from the number of
         * frames and frame rate.  But for the purpose of this program,
         * it's not necessary.
         */
        public Time getDuration()
        {
            return DURATION_UNKNOWN;
        }

        public Object[] getControls()
        {
            return new Object[0];
        }

        public Object getControl(String type)
        {
            return null;
        }
    }


    /**
     * The source stream to go along with ImageDataSource.
     */
    class ImageSourceStream implements PullBufferStream
    {
        Vector images;
        int width, height;
        VideoFormat format;
        StsProgressPanel panel = null;
        double progressStep = 1.0;

        int nextImage = 0;	// index of the next image to be read.
        boolean ended = false;

        public ImageSourceStream(int width, int height, int frameRate, Vector images, StsProgressPanel panel)
        {
            this.width = width;
            this.height = height;
            this.images = images;
            String imageFile = (String)images.elementAt(nextImage);
            String fmt = VideoFormat.JPEG;
            if(imageFile.contains("bmp"))
            	fmt = VideoFormat.RGB;
            	
            this.panel = panel;
            progressStep = 1.0/(images.size()+1);
            format = new VideoFormat(fmt, new Dimension(width, height),
                                     Format.NOT_SPECIFIED,
                                     Format.byteArray,
                                     (float)frameRate);
        }

        /**
         * We should never need to block assuming data are read from files.
         */
        public boolean willReadBlock()
        {
            return false;
        }

        /**
         * This is called from the Processor to read a frame worth
         * of video data.
         */
        public void read(Buffer buf) throws IOException
        {
            // Check if we've finished all the frames.
            if (nextImage >= images.size())
            {
                // We are done.  Set EndOfMedia.
 //               panel.appendLine("Done reading all images.");
                buf.setEOM(true);
                buf.setOffset(0);
                buf.setLength(0);
                ended = true;
                return;
            }

            String imageFile = (String)images.elementAt(nextImage);
            nextImage++;

            panel.appendLine("Reading image file: " + imageFile);
            panel.setValue(nextImage * progressStep);

//            panel.appendLine("Reading image file: " + imageFile);
//            panel.setProgressBarValue(nextImage/images.size());

            // Open a random access file for the next image.
            RandomAccessFile raFile;
            raFile = new RandomAccessFile(imageFile, "r");

            byte data[] = null;

            // Check the input buffer type & size.

            if (buf.getData() instanceof byte[])
                data = (byte[])buf.getData();

            // Check to see the given buffer is big enough for the frame.
            if (data == null || data.length < raFile.length())
            {
                data = new byte[(int)raFile.length()];
                buf.setData(data);
            }

            // Read the entire JPEG image from the file.
            raFile.readFully(data, 0, (int)raFile.length());

            buf.setOffset(0);
            buf.setLength((int)raFile.length());
            buf.setFormat(format);
            buf.setFlags(buf.getFlags() | buf.FLAG_KEY_FRAME);

            // Close the random access file.
            raFile.close();
        }

        /**
         * Return the format of each video frame.  That will be JPEG.
         */
        public Format getFormat()
        {
            return format;
        }

        public ContentDescriptor getContentDescriptor()
        {
            return new ContentDescriptor(ContentDescriptor.RAW);
        }

        public long getContentLength()
        {
            return 0;
        }

        public boolean endOfStream()
        {
            return ended;
        }

        public Object[] getControls()
        {
            return new Object[0];
        }

        public Object getControl(String type)
        {
            return null;
        }
    }

    /*
    public class JMFTest extends JFrame
    {
        Player _player;
        JMFTest()
        {
            addWindowListener( new WindowAdapter()
            {
                public void windowClosing( WindowEvent e )
                {
                    _player.stop();
                    _player.deallocate();
                    _player.close();
                    System.exit( 0 );
                }
            });
            setExtent( 0, 0, 320, 260 );
            JPanel panel = (JPanel)getContentPane();
            panel.setLayout( new BorderLayout() );
            String mediaFile = "vfw://1";
            try
            {
                MediaLocator mlr = new MediaLocator( mediaFile );
                _player = Manager.createRealizedPlayer( mlr );
                if (_player.getVisualComponent() != null)
                panel.add("Center", _player.getVisualComponent());
                if (_player.getControlPanelComponent() != null)
                panel.add("South", _player.getControlPanelComponent());
            }

            catch (Exception e)
            {
                System.err.println( "Got exception " + e );
            }
        }

        public static void main(String[] args)
        {
            JMFTest jmfTest = new JMFTest();
            jmfTest.show();
        }
    }
   */

    /*
    public class PlayerApplet extends Applet
    {
      Player player = null;

      public void init()
      {
        setLayout( new BorderLayout() );
        String mediaFile = getParameter( "FILE" );
        try
        {
          URL mediaURL = new URL( getDocumentBase(), mediaFile );
          player = Manager.createRealizedPlayer( mediaURL );
          if (player.getVisualComponent() != null)
            add("Center", player.getVisualComponent());

          if (player.getControlPanelComponent() != null)
            add("South", player.getControlPanelComponent());

        }
        catch (Exception e)
        {
          System.err.println( "Got exception " + e );
        }
      }

      public void start()
      {
        player.start();
      }

      public void stop()
      {
        player.stop();
        player.deallocate();
      }

      public void destroy()
      {
        player.close();
      }
    }
    */
}
