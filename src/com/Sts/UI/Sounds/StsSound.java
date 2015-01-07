package com.Sts.UI.Sounds;

import com.Sts.Utilities.*;
import sun.audio.*;

public class StsSound
{
    /** Available Sounds */
	public static String CAMERA_CLICK = "click.wav";
	public static String CAR_ALARM = "car_alarm.wav";    
    public static String MISSLEBEEP = "misslebeep.wav";
	public static String PING = "ping.wav";
	public static String RADARPING = "radarping.wav";
	public static String SINGLEGUNSHOT = "singlegunshot.wav";
	public static String THUNDER = "thunderrumble.wav";
	public static String SUBALERT = "submarinealert.wav";
	public static String IMPLOSION = "implosion2.wav";
	public static String HOMEALARM = "homealarm.wav";
	public static String M50 = "m50.wav";
	public static String GUNSHOT = "gunshots3.wav";
    public static String BOING = "boing.wav";
	public static String BEEP1 = "beep1.wav";
	public static String BEEP2 = "beep2.wav";
	public static String BEEP3 = "beep3.wav";
	public static String BEEP4 = "beep4.wav";
	public static String CAMERA = "camera.wav";	
	public static String CLICK = "click.wav";
	public static String CLANG = "clang.wav";
	public static String BUZZ = "buzz.wav";
	public static String KNOCK = "knock.wav";	
	public static String MINI2 = "mini2.wav";
	public static String SELECT1 = "select1.wav";
	public static String SELECT6 = "select6.wav";


    /** List of available sounds */
	public static String[] sounds = {CAMERA_CLICK, CAR_ALARM, MISSLEBEEP, PING, RADARPING,
            SINGLEGUNSHOT, THUNDER, SUBALERT, IMPLOSION, HOMEALARM,
            M50, GUNSHOT, BOING, BEEP1,
		BEEP2, BEEP3, BEEP4, CAMERA, CLICK, CLANG, BUZZ, KNOCK, MINI2, SELECT1, SELECT6};

    /**
     * Play a wave file
     * @param name - the name of the sound to play. must be a .wav file.
     */
	static public void play(String name)
	{
        play(StsSound.class, name);
    }
	static public void play(Class c, String name)
	{
		try
		{
			//java.net.URL url = new URL("Sounds/" + name);
            java.net.URL url = c.getResource(name);
			if (url == null) return;
			AudioStream as = new AudioStream (url.openStream());
			AudioPlayer.player.start(as);
//			JApplet.newAudioClip(url).play();
		}
		catch (Exception e)
		{
			StsException.outputException("Error playing audio file", e, StsException.WARNING);
        }
	}

    /**
     * Main to test StsSound object
     * @param args
     */
    public static void main(String[] args)
    {
		StsSound.play(BOING);
		StsSound.play(CAMERA_CLICK);
    }
}
