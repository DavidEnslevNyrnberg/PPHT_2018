package com.example.erik.acutestressapp;


import android.content.Intent;
import android.graphics.Color;
import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;
import android.view.View;
import android.widget.Button;
import android.widget.TextView;

import com.jjoe64.graphview.DefaultLabelFormatter;
import com.jjoe64.graphview.GraphView;
import com.jjoe64.graphview.helper.DateAsXAxisLabelFormatter;
import com.jjoe64.graphview.series.DataPoint;
import com.jjoe64.graphview.series.LineGraphSeries;

import com.opencsv.CSVReader;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.stat.StatUtils;

import java.io.IOException;
import java.io.InputStreamReader;
import java.text.SimpleDateFormat;
import java.util.Date;



public class MainActivity extends AppCompatActivity /*implements EmpaDataDelegate, EmpaStatusDelegate*/ {


	// For cleaner code. If one wants to connect to the E4 one must insert the content inside of
	// this function HERE
	//empatica_start()

	private TextView statusLabel;

	int lastX = 300;
	int lastXpeak = 297;

	double eda_fs;
	double bvp_fs;
	private double[] edaData;
	private double[] bvpData;
	private double[] bvpHRV;
	private double[] hrv_feature;
	java.util.Date[] time;
	private int[] bvp_peaks;
	private double[] eda_peaks;
	boolean pause = false;
	boolean start = false;
	boolean showeda = false;


	// Datapoint variables for the graph
	private LineGraphSeries<DataPoint> series_EDA;
	private LineGraphSeries<DataPoint> series_BVP;
	private LineGraphSeries<DataPoint> series_peak;
	SimpleDateFormat mmss = new SimpleDateFormat("hh:mm:ss");


	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_main);


		// Initialize vars that reference UI components
		statusLabel = findViewById(R.id.status);

		// Initializing the E4 device manager
		//initEmpaticaDeviceManager();




		// initial signal processing
		edaData = getEDA();
		bvpData = getBVP();

		bvp_peaks = findpeaks(bvpData, 0, 10);
		//eda_peaks = findpeaks(edaData,0.1f,6);
		bvpHRV = getHRV(bvp_peaks, bvpData[1], 0);
		hrv_feature = hrv_features(bvpHRV, 0, 15);
		time = time_to_timestamp(bvpData);

		// Runs the function that initializes the graphview and adds the EDA and BVP series
		initgraph();
	}



	private void addEDAEntry() {
	/*=======================================================================================
	EDA and peak data is appended to as a new datapoint on the graph. It's using a integer
	value LastX that just iterates from every fatapoint added. When adding a new datapoint
	it will scroll to the end, and will delete old datapoints from the graph that is
	longer than 10 000 points old.
	=======================================================================================*/

		series_EDA.appendData(new DataPoint(time[-300+lastX++], edaData[lastX]), true, 10000);
	}

	private void addBVPEntry() {
	/*=======================================================================================
	BVP and peak data is appended to as a new datapoint on the graph. It's using a integer
	value LastX that just iterates from every datapoint added. When adding a new datapoint
	it will scroll to the end, and will delete old datapoints from the graph that is
	longer than 10 000 points old. The peak series has an offset in the code, so thats
	taken care of setting an integer value with an offset of 3 samples.
	=======================================================================================*/

		series_BVP.appendData(new DataPoint(time[-300+lastX++], bvpData[lastX]), true,
		10000);
		series_peak.appendData(new DataPoint(time[-300+lastX], bvp_peaks[lastXpeak++]),
		true, 10000);


		/* Proof of concept testing of the Alert view update is working
		if(lastX>4000){
			updateTextView("Stressed",2);


		}else{
			updateTextView("Not stressed",1);
		}
		*/
	}

	private double[] getEDA() {
	/*=======================================================================================
	This function is reading the eda.csv file from the assets folder and appending the data
	to an array, is sent through a low pass filter with a cut-off frequency of 15, and put in
	to a new array with the correct size corresponding to the filtered signal

	OUTPUT: double array with a filtered eda signal
	=======================================================================================*/
		double[] res2 = new double[10000];

		int i = 0;

		try {

			// Reading the CSV file
			CSVReader reader = new CSVReader(new InputStreamReader(getAssets().open("eda.csv")));
			String[] next;
			while ((next = reader.readNext()) != null) {
				res2[i] = Double.parseDouble(next[0]);
				i += 1;
			}
			eda_fs = res2[1];


		} catch (IOException e) {
			e.printStackTrace();
		}
		double[] res = new double[i];
		double[] filter = new double[10000];
		System.arraycopy(res2, 2, filter, 2, res2.length - 2);
		filter = LowPassFilter(filter, 15, eda_fs);
		System.arraycopy(res2, 0, res, 0, 2);
		System.arraycopy(filter, 2, res, 2, res.length - 2);

		return res;
	}

	private double[] getBVP() {
	/*=======================================================================================
	This function is reading the bvp.csv file from the assets folder and appending the data
	to an array, is sent through a low pass filter with a cut-off frequency of 15, and put in
	to a new array with the correct size corresponding to the filtered signal

	OUTPUT: double array with a filtered bvp signal
	=======================================================================================*/
		double[] res2 = new double[131072];
		int j = 0;

		try {

			// CSV read
			CSVReader reader = new CSVReader(new InputStreamReader(getAssets().open("bvp.csv")));
			String[] next;
			while ((next = reader.readNext()) != null) {
				res2[j] = Double.parseDouble(next[0]);
				j += 1;
			}
			bvp_fs = res2[1];

		} catch (IOException e) {
			e.printStackTrace();
		}
		double[] res = new double[j];


		// LP filtering
		// TODO: Fix the low pass filtering of the signal


		double[] filter = new double[131072];
		System.arraycopy(res2, 2, filter, 2, res2.length - 2);
		filter = LowPassFilter(filter, 15, bvp_fs);
		System.arraycopy(res2, 0, res, 0, 2);
		System.arraycopy(filter, 2, res, 2, res.length - 2);
		return res;
	}

	// Update text in Stress status
	public void updateTextView(String text, int color) {
	/*=======================================================================================
	Updates the text and color of the stress status textview

	INPUTS: new string to show in the Textview, color of the background of the textview
	=======================================================================================*/
		TextView textView = statusLabel;
		textView.setText(text);
		switch (color) {
			case 1:
				textView.setBackgroundColor(0xFF19B796);
				break;
			case 2:
				textView.setBackgroundColor(Color.RED);

		}

	}

	public Date[] time_to_timestamp(double[] data) {
	/*=======================================================================================
	Takes the timestamp from the data array and converts it to Java Date format to be used in
	the x-axis on the graph

	INPUT: Signal data to get the time for the x-axis from
	OUTPUT: Array of the timestamp corresponding to the sample from the input
	=======================================================================================*/

		int start_stamp = (int) data[0];

		Date[] dates = new Date[data.length];


		for (int i = 0; i < data.length; i++) {
		// Makes the unix timestamp and adds on the time since the first sample (in milliseconds).
			Long time_conv = ((Long.valueOf(start_stamp)*1000)+(i*16));
				Date time = new java.util.Date(time_conv);
				dates[i] = time;

		}
		return dates;
	}

	public void alarmHistory(View v) {
		// Starts the activity to view the alarm history. Uses an on click listener in the xml file.
		Intent i = new Intent(getApplicationContext(), alarmHistory.class);
		//StartActivity with the intent i
		startActivity(i);

	}

	public void settings(View v) {
	/*=======================================================================================
	Starts the activity to view the settings. Uses an on click listener in the xml file.
	It will also return with the updated settings values to be applied.
	=======================================================================================*/

		Intent i = new Intent(getApplicationContext(), settings.class);
		//Start the activity with the intent i and requestcode 1
		startActivityForResult(i, 1);

	}

	// Not implemented
	@Override
	protected void onActivityResult(int requestCode, int resultCode, Intent data) {
	/*=======================================================================================
	When the settings activity ends, this function will run and handle the return data from
	the changes in the settings activity.
	=======================================================================================*/
		super.onActivityResult(requestCode, resultCode, data);
		if (requestCode == 1 && resultCode == RESULT_OK) {
			showeda = getIntent().getExtras().getBoolean("data", false);

		}
	}

	public void startstop(View v) {
	/*=======================================================================================
	On click listener for the start stop button. When clicked it starts a new thread which
	appends BVP or EDA data to the graph, depending on what is selected to be viewed in the graph
	. First time pressed it will start the thread. If pressed more times it will toggle between
	pausing the appending of new data, or continuing appending data points to the graph
	=======================================================================================*/

		Button startstop_button = (Button) findViewById(R.id.startStop);
		if (!start) {
			start = true;
			startstop_button.setText("Pause");

			new Thread(new Runnable() {

				@Override
				public void run() {
					// One problem that will occur here is that the length of the EDA is shorter
					// than the BVP array, and will cause an overflow.
					for (int i = 0; i < bvpData.length; i++) {
						runOnUiThread(new Runnable() {

							@Override
							public void run() {

								if (!pause) {
									if (showeda) {
										addEDAEntry();
									} else {
										addBVPEntry();
									}
								}
							}
						});


						// sleep function that corresponds to the sample frequency of the signal
						// data
						try {
							if (showeda) {
								Thread.sleep(250);
							} else {
								Thread.sleep(15);
							}
						} catch (InterruptedException e) {
							// manage error ...
						}
					}
					start = false;
				}
			}).start();

		} else {
			if (!pause) {
				pause = true;
				startstop_button.setText("Resume");

			} else {
				pause = false;
				startstop_button.setText("Pause");
			}
		}
	}

	public double[] LowPassFilter(double[] data, double lowPass, double frequency) {
		/*=======================================================================================
		source: https://stackoverflow.com/questions/4026648/how-to-implement-low-pass-filter
		data: input data, must be spaced equally in time.
		lowPass: The cutoff frequency at which
		frequency: The frequency of the input data.

		The apache FFT accepts arrays that are powers of 2. and therefore the array must be
		padded with zeros to get the correct length.
		=======================================================================================*/
		int minPowerOf2 = 1;
		while (minPowerOf2 < data.length)
			minPowerOf2 = 2 * minPowerOf2;

		//pad with zeros
		double[] padded = new double[minPowerOf2];
		for (int i = 0; i < data.length; i++)
			padded[i] = data[i];


		FastFourierTransformer transformer = new FastFourierTransformer(DftNormalization.STANDARD);
		Complex[] fourierTransform = transformer.transform(padded, TransformType.FORWARD);

		//build the frequency domain array
		double[] frequencyDomain = new double[fourierTransform.length];
		for (int i = 0; i < frequencyDomain.length; i++)
			frequencyDomain[i] = frequency * i / (double) fourierTransform.length;

		//build the classifier array, 2s are kept and 0s do not pass the filter
		double[] keepPoints = new double[frequencyDomain.length];
		keepPoints[0] = 1;
		for (int i = 1; i < frequencyDomain.length; i++) {
			if (frequencyDomain[i] < lowPass)
				keepPoints[i] = 2;
			else
				keepPoints[i] = 0;
		}

		//filter the fft
		for (int i = 0; i < fourierTransform.length; i++)
			fourierTransform[i] = fourierTransform[i].multiply((double) keepPoints[i]);

		//invert back to time domain
		Complex[] reverseFourier = transformer.transform(fourierTransform, TransformType.INVERSE);

		//get the real part of the reverse
		double[] result = new double[data.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = reverseFourier[i].getReal();
		}

		return result;
	}

	public int[] findpeaks(double[] data, float threshold, int window_size) {
		/*=======================================================================================
		This function takes the filtered data signal and checks for the local minima in small
		windows set in the input window. It also uses a manual threshold value to look for the
		local minima in. When it has iterated through the whole signal it returns the peak
		locations.

		INPUT: Filtered signal, threshold value, windows size
		OUTPUT: Signal with same length as input signal with the location of the peaks
		=======================================================================================*/

		int[] peaks = new int[data.length];

		int i = 300;
		int increment_size = window_size / 2;


		while (i < data.length - (data.length / 10)) {
			float avg_front = 0;
			float avg_mid = 0;
			float avg_back = 0;

			for (int k = 0; k < window_size; k++) {
				// Averaging values
				avg_mid += data[i + k] / window_size;
				avg_front += data[i + k - 5] / window_size;
				avg_back += data[i + k + 5] / window_size;
			}
			if (avg_back > avg_mid && avg_mid < avg_front && avg_mid < threshold) {
				peaks[i] = 50;
			} else {
				peaks[i] = 0;
			}
			i += increment_size;
		}

		return peaks;
	}

	public double[] getHRV(int[] data_peaks, double fs, int window_iteration) {
		// This function takes the array of the peak location and calculates the distance between
		// the peaks in the BVP signal in a window of 15 seconds. The window_iteration is used to
		// move though the data peaks signal in chunks of 15 seconds. The HRV is then
		// interpolated before returned
		//
		// INPUT: peaks location, sample frequency, what 15 second chunk of the signal is calculated
		// OUTPUT: interpolated HRV signal
		int window_pos = window_iteration * 15;
		int window_size = 15 * (int) fs;

		double[] init_HRV = new double[window_size];

		int peak_pos = 0;
		int hrv_iteration = 0;

		for (int i = 0; i < init_HRV.length; i++) {

			if (data_peaks[i] != 0) {

				peak_pos = i - peak_pos;
				init_HRV[hrv_iteration] = (double) peak_pos;
				peak_pos = i;
				hrv_iteration += 1;
			}
		}

		double[] HRV = new double[hrv_iteration * 8];

		int inter_iteration = 0;
		for (int i = 0; i < hrv_iteration - 5; i++) {
			double inter_step = (init_HRV[i + 1] - init_HRV[i]) / 8;
			for (int j = 1; j < 9; j++) {
				if (inter_step != 0f) {
					HRV[inter_iteration + j] = init_HRV[i] + (j * inter_step);
				} else {
					HRV[inter_iteration + j] = init_HRV[i];
				}
			}
			if (i == 1140) {
				inter_iteration += 1;
				inter_iteration -= 1;
			}
			inter_iteration += 8;
		}
		return HRV;
	}


	// Not implemented
	public double[] hrv_features(double[] HRV_signal, int window_iteration, int window_size) {
		/*=======================================================================================
		This function extract the features of the HRV signal described in the project journal.
		The window_iteration is used to move though the HRV signal in chunks of
		15 seconds. The window size can be manually set, but since the HRV signal is set to 15
		seconds the HRV is windowed to 15 seconds that is the max size of the window.

		INPUT: HRV signal,  what 15 second chunk of the signal is calculated, size of the window
		=======================================================================================*/
		/*=======================================================================================

		=======================================================================================*/

		int fs = 8;
		window_size = window_size * 8;
		FastFourierTransformer FFT = new FastFourierTransformer(DftNormalization.STANDARD);

		int iteration_multiplier = (window_iteration * 16); // -1 to get the signal to a power of 2
		double[] hrv_windowed = new double[128];
		double[] hrv_features = new double[11];
		System.arraycopy(HRV_signal, iteration_multiplier, hrv_windowed, 0, HRV_signal.length);

		double mean = StatUtils.mean(hrv_windowed);
		double sd_nn = get_std(hrv_windowed);
		double rms_nn = rootMeanSquare(hrv_windowed);
		double sdsd = get_std(get_diff(hrv_windowed));
		double nn50 = get_nn50(get_diff(hrv_windowed));
		double pnn50 = (nn50 / hrv_windowed.length * 100);
		// FREQ
		Complex[] HRV_fft = FFT.transform(hrv_windowed, TransformType.FORWARD);


		double HRV_power = HRV_power(HRV_fft);
		// Not implemented
		double area_vlf;
		double area_lf_norm;
		double area_hf_norm;
		double lf_hf_ratio;

		hrv_features[0] = mean;
		hrv_features[1] = sd_nn;
		hrv_features[2] = rms_nn;
		hrv_features[3] = sdsd;
		hrv_features[4] = nn50;
		hrv_features[5] = pnn50;
		hrv_features[6] = HRV_power;
		/* Not calculated yet
		hrv_features[7] = area_vlf;
		hrv_features[8] = area_lf_norm;
		hrv_features[9] = area_hf_norm;
		hrv_features[10] = lf_hf_ratio;
		*/
		return hrv_features;
	}


	public void SVM(double[] features){
		//TODO: IMPLEMENT THE SVM MODEL
	}

	public double HRV_power(Complex[] signal) {
		// Calculates the power off the fft signal
		//
		// INPUT: a fft signal
		// OUTPUT: power of the input signal
		double power = 0;
		for (int i = 0; i < signal.length; i++) {
			power += signal[i].abs();
		}
		power = Math.pow(power, 2);

		power = power / signal.length;

		return power;

	}

	// Not implemented
	public double[] evaluate_eda(double[] eda_signal, int window_iteration) {
		//TODO: Features in 15s windows 10
		double[] foo = new double[15];

		return foo;
	}



	// features calculation functions:
	// The mathematically formulas is given in the journal paper

	public static double rootMeanSquare(double[] nums) {
	/*Calculates the root mean square of the input array.
	*
	* INPUT: Array to */
		double sum = 0.0;
		for (double num : nums)
			sum += num * num;
		return Math.sqrt(sum / nums.length);
	}

	public double get_std(double[] hrv) {
		DescriptiveStatistics stats = new DescriptiveStatistics();
		for (int i = 0; i < hrv.length; i++) {
			stats.addValue(hrv[i]);
		}
		double res = stats.getStandardDeviation();
		return res;
	}


	public double[] get_diff(double[] signal) {
		double[] diff_signal = new double[signal.length];

		for (int i = 0; i < diff_signal.length - 1; i++) {
			diff_signal[i] = signal[i] - signal[i + 1];
		}
		return diff_signal;
	}

	public double get_nn50(double[] signal) {
		double nn50 = 0;
		int msInt = 50 / 1000;
		double[] over_msInt = new double[signal.length];

		for (int i = 0; i < signal.length; i++) {
			if (signal[i] != 0 && signal[i] > msInt) {
				nn50 += over_msInt[i];
			}
		}
		return nn50;
	}

	public void initgraph() {
		// This function is the initialization of the graph and the data series. It enables
		// scrolling and zooming of the graph, as well as setting the min and max values of the
		// graph window in both the x and y axis.

		GraphView graph = (GraphView) findViewById(R.id.graph);
		// Settings for the interaction with the graph
		// enables horizontal zooming and scrolling
		graph.getViewport().setScalable(true);
		// enables vertical zooming and scrolling
		graph.getViewport().setScalableY(true);
		// set manual Y bounds
		graph.getViewport().setYAxisBoundsManual(true);
		if (showeda) {
			graph.getViewport().setMinY(-1);
			graph.getViewport().setMaxY(1);
		} else {
			graph.getViewport().setMinY(-150);
			graph.getViewport().setMaxY(150);
		}
		graph.getViewport().setMinX(time[0].getTime());
		graph.getViewport().setMaxX(time[300].getTime());


		// Changes the x-axis of the graph window to display the time in HH:MM:SS
		graph.getGridLabelRenderer().setLabelFormatter(new DefaultLabelFormatter()
		{
			@Override
			public String formatLabel(double value, boolean isValueX){
				if(isValueX)
				{
					return mmss.format(new Date((long) value));
				}else {
				return super.formatLabel(value, isValueX);
				}
			}
		});

		graph.getGridLabelRenderer().setNumHorizontalLabels(4);
		graph.getGridLabelRenderer().setHumanRounding(false);


		// Initialization of the graph series
		// TODO: Fix the bug where the showing graph has to be declared first
		series_BVP = new LineGraphSeries<>();
		graph.addSeries(series_BVP);
		series_peak = new LineGraphSeries<>();
		graph.addSeries(series_peak);
		series_EDA = new LineGraphSeries<>();
		graph.addSeries(series_EDA);


		// styling BVP graph
		series_BVP.setColor(Color.argb(255, 255, 0, 0));   // Setting the color
		series_BVP.setDrawDataPoints(false);                                // Disable dots in the graph drawing
		series_BVP.setThickness(3);                                         // Setting the thickness of the graph series
		series_BVP.setDrawAsPath(true);                                     // Drawing as a path

		// styling EDA series the same way as the BVP series
		series_EDA.setColor(Color.argb(255, 0, 255, 0));
		series_EDA.setDrawDataPoints(false);
		series_EDA.setThickness(3);
		// as we use dates as labels, the human rounding to nice readable numbers
		// is not necessary

		// styling peak series the same way as the BVP series
		series_peak.setColor(Color.argb(255, 0, 200, 255));
		series_peak.setDrawDataPoints(false);
		series_peak.setThickness(3);
	}

}

// Functions that would be used to use the E4 directly in live mode.
// Since live mode is not used in the application at this moment it is commented out and
// moved to under the program to get cleaner code


/*
	public void empatica_start(){
		// For E4 set up
		private static final int REQUEST_ENABLE_BT = 1;
		private static final int REQUEST_PERMISSION_ACCESS_COARSE_LOCATION = 1;
		private static final long STREAMING_TIME = 100000; // Stops streaming 100 000 seconds after
		// connection
		private static final String EMPATICA_API_KEY = "2baa364dceb843299b942b41a276740b";
		private EmpaDeviceManager deviceManager = null;
	}


	@Override
	public void onRequestPermissionsResult(int requestCode, @NonNull String[] permissions, @NonNull int[] grantResults) {
		switch (requestCode) {
			case REQUEST_PERMISSION_ACCESS_COARSE_LOCATION:
				// If request is cancelled, the result arrays are empty.
				if (grantResults.length > 0 && grantResults[0] == PackageManager.PERMISSION_GRANTED) {
					// Permission was granted, yay!
					initEmpaticaDeviceManager();
				} else {
					// Permission denied, boo!
					final boolean needRationale = ActivityCompat.shouldShowRequestPermissionRationale(this, Manifest.permission.ACCESS_COARSE_LOCATION);
					new AlertDialog.Builder(this)
							.setTitle("Permission required")
							.setMessage("Without this permission bluetooth low energy devices cannot be found, allow it in order to connect to the device.")
							.setPositiveButton("Retry", new DialogInterface.OnClickListener() {
								public void onClick(DialogInterface dialog, int which) {
									// try again
									if (needRationale) {
										// the "never ask again" flash is not set, try again with permission request
										initEmpaticaDeviceManager();
									} else {
										// the "never ask again" flag is set so the permission requests is disabled, try open app settings to enable the permission
										Intent intent = new Intent(Settings.ACTION_APPLICATION_DETAILS_SETTINGS);
										Uri uri = Uri.fromParts("package", getPackageName(), null);
										intent.setData(uri);
										startActivity(intent);
									}
								}
							})
							.setNegativeButton("Exit application", new DialogInterface.OnClickListener() {
								public void onClick(DialogInterface dialog, int which) {
									// without permission exit is the only way
									finish();
								}
							})
							.show();
				}

				break;
		}

	private void initEmpaticaDeviceManager() {
		// Android 6 (API level 23) now require ACCESS_COARSE_LOCATION permission to use BLE
		if (ContextCompat.checkSelfPermission(this, Manifest.permission.ACCESS_COARSE_LOCATION) != PackageManager.PERMISSION_GRANTED) {
			ActivityCompat.requestPermissions(this, new String[] { Manifest.permission.ACCESS_COARSE_LOCATION }, REQUEST_PERMISSION_ACCESS_COARSE_LOCATION);
		} else {
			// Create a new EmpaDeviceManager. MainActivity is both its data and status delegate.
			deviceManager = new EmpaDeviceManager(getApplicationContext(), this, this);

			if (TextUtils.isEmpty(EMPATICA_API_KEY)) {
				new AlertDialog.Builder(this)
						.setTitle("Warning")
						.setMessage("Please insert your API KEY")
						.setNegativeButton("Close", new DialogInterface.OnClickListener() {
							public void onClick(DialogInterface dialog, int which) {
								// without permission exit is the only way
								finish();
							}
						})
						.show();
				return;
			}
			// Initialize the Device Manager using your API key. You need to have Internet access at this point.
			deviceManager.authenticateWithAPIKey(EMPATICA_API_KEY);
		}
	}

	@Override
	protected void onPause() {
		super.onPause();
		if (deviceManager != null) {
			deviceManager.stopScanning();
		}
	}

	@Override
	protected void onDestroy() {
		super.onDestroy();
		if (deviceManager != null) {
			deviceManager.cleanUp();
		}
	}

	@Override
	public void didDiscoverDevice(BluetoothDevice bluetoothDevice, String deviceName, int rssi, boolean allowed) {
		// Check if the discovered device can be used with your API key. If allowed is always false,
		// the device is not linked with your API key. Please check your developer area at
		// https://www.empatica.com/connect/developer.php
		if (allowed) {
			// Stop scanning. The first allowed device will do.
			deviceManager.stopScanning();
			try {
				// Connect to the device
				deviceManager.connectDevice(bluetoothDevice);
				updateLabel(deviceNameLabel, "To: " + deviceName);
			} catch (ConnectionNotAllowedException e) {
				// This should happen only if you try to connect when allowed == false.
				Toast.makeText(MainActivity.this, "Sorry, you can't connect to this device", Toast.LENGTH_SHORT).show();
			}
		}
	}

	@Override
	public void didRequestEnableBluetooth() {
		// Request the user to enable Bluetooth
		Intent enableBtIntent = new Intent(BluetoothAdapter.ACTION_REQUEST_ENABLE);
		startActivityForResult(enableBtIntent, REQUEST_ENABLE_BT);
	}

	@Override
	protected void onActivityResult(int requestCode, int resultCode, Intent data) {
		// The user chose not to enable Bluetooth
		if (requestCode == REQUEST_ENABLE_BT && resultCode == Activity.RESULT_CANCELED) {
			// You should deal with this
			return;
		}
		super.onActivityResult(requestCode, resultCode, data);
	}

	@Override
	public void didUpdateSensorStatus(EmpaSensorStatus status, EmpaSensorType type) {
		// No need to implement this right now
	}

	@Override
	public void didUpdateStatus(EmpaStatus status) {
		// Update the UI
		updateLabel(statusLabel, status.name());

		// The device manager is ready for use
		if (status == EmpaStatus.READY) {
			updateLabel(statusLabel, status.name() + "");
			// Start scanning
			deviceManager.startScanning();
			// The device manager has established a connection
		} else if (status == EmpaStatus.CONNECTED) {
			// Stop streaming after STREAMING_TIME
			runOnUiThread(new Runnable() {
				@Override
				public void run() {
					dataCnt.setVisibility(View.VISIBLE);
					new Handler().postDelayed(new Runnable() {
						@Override
						public void run() {
							// Disconnect device
							deviceManager.disconnect();
						}
					}, STREAMING_TIME);
				}
			});
			// The device manager disconnected from a device
		} else if (status == EmpaStatus.DISCONNECTED) {
			updateLabel(deviceNameLabel, "");
		}
	}


	@Override
	public void didReceiveAcceleration(int x, int y, int z, double timestamp) {
		//updateLabel(accel_xLabel, "" + x);
		//updateLabel(accel_yLabel, "" + y);
		//updateLabel(accel_zLabel, "" + z);
	}

	@Override
	public void didReceiveBVP(float bvp, double timestamp) {
		//updateLabel(bvpLabel, "" + bvp);

	}

	@Override
	public void didReceiveBatteryLevel(float battery, double timestamp) {
		//updateLabel(batteryLabel, String.format("%.0f %%", battery * 100));
	}

	@Override
	public void didReceiveGSR(float gsr, double timestamp) {
		updateLabel(edaLabel, "" + gsr);
	}

	@Override
	public void didReceiveIBI(float ibi, double timestamp) {
		//updateLabel(ibiLabel, "" + ibi);
	}

	@Override
	public void didReceiveTemperature(float temp, double timestamp) {
		//updateLabel(temperatureLabel, "" + temp);
	}

	// Update a label with some text, making sure this is run in the UI thread
	private void updateLabel(final TextView label, final String text) {
		//runOnUiThread(new Runnable() {
			@Override
			public void run() {
				label.setText(text);
			}
		});
	}
}
}*/









