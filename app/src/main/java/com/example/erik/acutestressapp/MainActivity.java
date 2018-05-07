package com.example.erik.acutestressapp;


import android.content.Intent;
import android.graphics.Color;
import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;

import android.view.View;
import android.widget.Button;
import android.widget.TextView;

import com.jjoe64.graphview.GraphView;
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
import java.util.Date;



public class MainActivity extends AppCompatActivity /*implements EmpaDataDelegate, EmpaStatusDelegate*/ {

	// For E4 set up
	/*private static final int REQUEST_ENABLE_BT = 1;
	private static final int REQUEST_PERMISSION_ACCESS_COARSE_LOCATION = 1;

	private static final long STREAMING_TIME = 100000; // Stops streaming 100 000 seconds after
	// connection

	private static final String EMPATICA_API_KEY = "2baa364dceb843299b942b41a276740b";

	private EmpaDeviceManager deviceManager = null;*/




	//private TextView accel_xLabel;
	//private TextView accel_yLabel;
	//private TextView accel_zLabel;
	//private TextView bvpLabel;
	//private TextView edaLabel;
	//private TextView ibiLabel;
	//private TextView temperatureLabel;
	//private TextView batteryLabel;
	private TextView statusLabel;
	//private TextView deviceNameLabel;
	//private RelativeLayout dataCnt;



	int lastX = 300;

	double eda_fs;
	double bvp_fs;
	private Float[] edaData;
	private Float[] bvpData;
	private double[] bvpHRV;
	private double[] evaluated_hrv;
	java.util.Date time;
	private int[] bvp_peaks;
	private Float[] eda_peaks;
	boolean pause = false;
	boolean start = false;
	boolean showeda = true;


	// Datapoint variables for the graph
	private LineGraphSeries<DataPoint> series_EDA;
	private LineGraphSeries<DataPoint> series_BVP;


	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_main);


		// Initialize vars that reference UI components
		statusLabel =  findViewById(R.id.status);
		//dataCnt =  findViewById(R.id.dataArea);
		//accel_xLabel =  findViewById(R.id.accel_x);
		//accel_yLabel =  findViewById(R.id.accel_y);
		//accel_zLabel =  findViewById(R.id.accel_z);
		//bvpLabel =  findViewById(R.id.bvp);
		//edaLabel =  findViewById(R.id.eda);
		//ibiLabel =  findViewById(R.id.ibi);
		//temperatureLabel =  findViewById(R.id.temperature);
		//batteryLabel =  findViewById(R.id.battery);
		//deviceNameLabel =  findViewById(R.id.deviceName);
		GraphView graph = (GraphView) findViewById(R.id.graph);


		// Initializing the E4
		//initEmpaticaDeviceManager();

		// Settings for the interaction with the graph
		// enables horizontal zooming and scrolling
		graph.getViewport().setScalable(true);
		// enables vertical zooming and scrolling
		graph.getViewport().setScalableY(true);
		// set manual Y bounds
		graph.getViewport().setYAxisBoundsManual(true);
		if(showeda){
			graph.getViewport().setMinY(-1);
			graph.getViewport().setMaxY(1);
			graph.getGridLabelRenderer().setNumHorizontalLabels(3);
			graph.getGridLabelRenderer().setHumanRounding(true);
		}else{
			graph.getViewport().setMinY(-150);
			graph.getViewport().setMaxY(150);
			graph.getGridLabelRenderer().setNumHorizontalLabels(3);
			graph.getGridLabelRenderer().setHumanRounding(true);
		}
		graph.getViewport().setMinX(300);
		graph.getViewport().setMaxX(500);



		// Initialization of the graph series
		// TODO: Fix the bug where the showing graph has to be declared first
		series_EDA = new LineGraphSeries<>();
		graph.addSeries(series_EDA);
		series_BVP = new LineGraphSeries<>();
		graph.addSeries(series_BVP);


		// styling BVP graph
		// Setting the color
		series_BVP.setColor(Color.argb(255,255,0,0));
		// Disable dots in the graph drawing
		series_BVP.setDrawDataPoints(false);
		// Setting the thickness of the graph series
		series_BVP.setThickness(3);
		// Drawing as a path
		series_BVP.setDrawAsPath(true);
		// styling EDA series the same way as the BVP series
		series_EDA.setColor(Color.argb(255,0,255,0));
		series_EDA.setDrawDataPoints(false);
		series_EDA.setThickness(3);
		// as we use dates as labels, the human rounding to nice readable numbers
		// is not necessary



		// initial signal processing
		edaData = getEDA();
		bvpData = getBVP();

		bvp_peaks = findpeaks(bvpData,0, 10);
		//eda_peaks = findpeaks(edaData,0.1f,6);
		bvpHRV = getHRV(bvp_peaks, bvpData[1], 0);
		evaluated_hrv = evaluate_hrv(bvpHRV,0,16);

	}




	// Appending data to graph
	private void addEDAEntry() {
		// here, we choose to display max 1000 points on the viewport and we scroll to end

		series_EDA.appendData(new DataPoint(lastX++, edaData[lastX]), true, 1000);
	}
	private void addBVPEntry() {

		//here, we choose to display max 1000 points on the viewport and we scroll to end

		series_BVP.appendData(new DataPoint(lastX++, bvpData[lastX]), true, 10000);

		/* Proof of concept testing of the Alert view update is working
		if(lastX>4000){
			updateTextView("Stressed",2);


		}else{
			updateTextView("Not stressed",1);
		}
		*/
	}

	// Function to read the EDA file and do the pre-processing
	// TODO: Finish pre-processing
	private Float[] getEDA(){
		/*This function is reading the eda.csv file from the assets folder and appending the data
		 to an array for further processing.

		 OUTPUT: Float array with eda signal
		 */
		Float[] res2 = new Float[10000];

		int i = 0;

		try {

			// Reading the CSV file
			CSVReader reader = new CSVReader(new InputStreamReader(getAssets().open("eda.csv")));
			String[] next;
			while ((next = reader.readNext()) != null){
				res2[i] = Float.parseFloat(next[0]);
				i += 1;
			}
			eda_fs = res2[1];

		} catch (IOException e) {
			e.printStackTrace();
		}

		Float[] res = new Float[i];
			System.arraycopy(res2,0,res,0,res.length);
		// LP filtering
		// TODO: Fix the low pass filtering of the signal
		//res = fourierLowPassFilter(res, 15, eda_fs);
		return res;
	}
	// Function to read the BVP file and do the pre-processing
	// TODO: Finish pre-processing
	private Float[] getBVP(){
		/*This function is reading the bvp.csv file from the assets folder and appending the data
		 to an array for further processing.

		 OUTPUT: Float array with bvp signal
		 */
		Float[] res2 = new Float[100000];
		int j = 0;

		try {

			// CSV read
			CSVReader reader = new CSVReader(new InputStreamReader(getAssets().open("bvp.csv")));
			String[] next;
			while ((next = reader.readNext()) != null){
				res2[j] = Float.parseFloat(next[0]);
				j += 1;
			}
			bvp_fs = res2[1];

		} catch (IOException e) {
			e.printStackTrace();
		}
		Float[] res = new Float[j];
		System.arraycopy(res2,0,res,0,res.length);

		// LP filtering
		// TODO: Fix the low pass filtering of the signal
		//res2 = fourierLowPassFilter(res2, 15, bvp_fs);
		return res;
	}

	// Update text in Stress status
	public void updateTextView(String text, int color) {
	/* Updates the text and color of the stress status textview
	*
	* INPUTS: new string to show in the Textview, color of the background of the textview*/
		TextView textView = statusLabel;
		textView.setText(text);
		switch (color){
			case 1:
				textView.setBackgroundColor(0xFF19B796);
				break;
			case 2:
				textView.setBackgroundColor(Color.RED);

		}

	}
	// Not implemented
	public Date[] time_to_timestamp(Float[] data){
	/* Takes the timestamp from the data array and converts it to Java Date format
	*
	* INPUT: Signal data to get the time for the x-axis from*/

		int start_stamp = data[0].intValue();
		int fs = data[1].intValue();
		int timestamp = 0;

		Date[] dates = new Date[data.length/fs];


		for(int i=0; i<data.length-fs; i++){
			if((i%fs) == 0){

				Date time = new java.util.Date((long)(start_stamp+(i/fs))*1000);
				dates[i/fs] = time;
			}

		}
		return dates;
	}

	public void alarmHistory(View v){
		// Starts the activity to view the alarm history. Uses an on click listener in the xml file.
		Intent i = new Intent(getApplicationContext(),alarmHistory.class);
		//StartActivity with the intent i
		startActivity(i);

	}
	public void settings(View v){
		/*Starts the activity to view the settings. Uses an on click listener in the xml file.
		* It will also return with the updated settings values to be applied. */
		Intent i = new Intent(getApplicationContext(),settings.class);
		//Start the activity with the intent i and requestcode 1
		startActivityForResult(i,1);

	}

	// Receiving the result from the StartActivityForResult()
	@Override
	protected void onActivityResult(int requestCode, int resultCode, Intent data) {
	/* when the settings activity ends, this function will run and handle the return data from
	the changes in the settings activity.*/
		super.onActivityResult(requestCode, resultCode, data);
		if (requestCode == 1 && resultCode == RESULT_OK){
			 showeda = getIntent().getExtras().getBoolean("data",false);

		}
	}

	public void startstop(View v){
	/*On click listener for the start stop button. When clicked it starts a new thread which
	appends BVP or EDA data to the graph, depending on what */

		Button startstop_button = (Button)findViewById(R.id.startStop);
		if(!start){
			start = true;
			startstop_button.setText("Pause");

			new Thread(new Runnable() {

				@Override
				public void run() {
					// we add the entries to the graph

					for (int i = 0; i < bvpData.length - 50; i++) { // TODO: fix length of i
						runOnUiThread(new Runnable() {

							@Override
							public void run() {

								if(!pause){
									if (showeda) {
										addEDAEntry();
									}else{
										addBVPEntry();
									}
								}



							}

						});


						// sleep to slow down the add of entries
						try {
							if(showeda){
								Thread.sleep(250);
							}else{
								Thread.sleep(16);
							}

						} catch (InterruptedException e) {
							// manage error ...
						}
					}
					start = false;

				}
			}).start();
		}else{
			if(!pause){
				pause = true;
				startstop_button.setText("Resume");

			}else{
				pause = false;
				startstop_button.setText("Pause");
			}
		}


	}
	// Not implemented
	public double[] fourierLowPassFilter(double[] data, double lowPass, double frequency){
		//data: input data, must be spaced equally in time.
		//lowPass: The cutoff frequency at which
		//frequency: The frequency of the input data.

		//The apache Fft (Fast Fourier Transform) accepts arrays that are powers of 2.
		int minPowerOf2 = 1;
		while(minPowerOf2 < data.length)
			minPowerOf2 = 2 * minPowerOf2;

		//pad with zeros
		double[] padded = new double[minPowerOf2];
		for(int i = 0; i < data.length; i++)
			padded[i] = data[i];


		FastFourierTransformer transformer = new FastFourierTransformer(DftNormalization.STANDARD);
		Complex[] fourierTransform = transformer.transform(padded, TransformType.FORWARD);

		//build the frequency domain array
		double[] frequencyDomain = new double[fourierTransform.length];
		for(int i = 0; i < frequencyDomain.length; i++)
			frequencyDomain[i] = frequency * i / (double)fourierTransform.length;

		//build the classifier array, 2s are kept and 0s do not pass the filter
		double[] keepPoints = new double[frequencyDomain.length];
		keepPoints[0] = 1;
		for(int i = 1; i < frequencyDomain.length; i++){
			if(frequencyDomain[i] < lowPass)
				keepPoints[i] = 2;
			else
				keepPoints[i] = 0;
		}

		//filter the fft
		for(int i = 0; i < fourierTransform.length; i++)
			fourierTransform[i] = fourierTransform[i].multiply((double)keepPoints[i]);

		//invert back to time domain
		Complex[] reverseFourier = transformer.transform(fourierTransform, TransformType.INVERSE);

		//get the real part of the reverse
		double[] result = new double[data.length];
		for(int i = 0; i< result.length; i++){
			result[i] = reverseFourier[i].getReal();
		}

		return result;
	}

	public int[] findpeaks(Float[] data, float threshold, int window_size){
		int[] peaks = new int[data.length];

		int i = 300;
		int increment_size = window_size/2;


		while(i < data.length - (data.length/10)){
			float avg_front = 0;
			float avg_mid = 0;
			float avg_back = 0;

			for(int k=0;k<window_size;k++){
				// Averaging values
				avg_mid += data[i+k]/window_size;
				avg_front += data[i+k-5]/window_size;
				avg_back += data[i+k+5]/window_size;
			}
			if( avg_back > avg_mid && avg_mid < avg_front && avg_mid < threshold){
				peaks[i] =  50;
			}
			else{
				peaks[i] = 0;
			}
			i += increment_size;
		}

		return peaks;
	}

	public double[] getHRV(int[] data_peaks, float fs, int window_iteration){
		int window_pos = window_iteration*16;
		int window_size = 16*(int)fs;

		double[] init_HRV = new double[window_size];


		int peak_pos = 0;
		int hrv_iteration = 0;


		for(int i = 0; i < init_HRV.length; i++){

			if(data_peaks[i] != 0){

				peak_pos = i - peak_pos;
				init_HRV[hrv_iteration] = (double) peak_pos;
				peak_pos = i;
				hrv_iteration +=1;
			}


		}

		double[] HRV = new double[hrv_iteration*8];

		int inter_iteration = 0;
		for(int i=0;i<hrv_iteration-5;i++){
			double inter_step = (init_HRV[i+1]-init_HRV[i])/8;
			for (int j=1;j<9;j++){
				if (inter_step != 0f){
					HRV[inter_iteration+j] = init_HRV[i]+(j*inter_step);
				}else{
					HRV[inter_iteration+j] = init_HRV[i];
				}

			}
			if(i==1140){
				inter_iteration +=1;
				inter_iteration -=1;
			}
			inter_iteration +=8;
		}
		return HRV;
	}


	// Not implemented
	public double[] evaluate_hrv(double[] HRV_signal, int window_iteration, int window_size){

		int fs = 8;
		window_size = window_size*8;
		FastFourierTransformer FFT = new FastFourierTransformer(DftNormalization.STANDARD);

		int iteration_multiplier = (window_iteration*16); // -1 to get the signal to a power of 2
		double[] hrv_windowed = new double[128];
		double[] hrv_features = new double[11];
		System.arraycopy(HRV_signal,iteration_multiplier,hrv_windowed,0,HRV_signal.length);

		double mean = StatUtils.mean(hrv_windowed);
		double sd_nn = get_std(hrv_windowed);
		double rms_nn = rootMeanSquare(hrv_windowed);
		double sdsd = get_std(get_diff(hrv_windowed));
		double nn50 = get_nn50(get_diff(hrv_windowed));
		double pnn50 = (nn50/hrv_windowed.length*100);
		// FREQ
		Complex[] HRV_fft = FFT.transform(hrv_windowed, TransformType.FORWARD);


		double HRV_power = HRV_power(HRV_fft);
		double area_vlf;
		double area_lf_norm;
		double area_hf_norm;
		double lf_hf_ratio;

		return hrv_features;
	}

	public double HRV_power(Complex[] signal){
		double power = 0;
		for (int i=0;i<signal.length;i++){
			power += signal[i].abs();
		}
		power = Math.pow(power,2);

		power = power/signal.length;

		return power;

	}
	// Not implemented
	public double[] evaluate_eda(double[] eda_signal, int window_iteration){
		double[] foo = new double[1];
		int fs = 4;
		int window_size  = 15;
		//TODO: Features in 15s windows 10
		//

		return foo;
	}

	//TODO: SVM mocdel that takes in evaluation of HRV and EDA

	// features calculation functions:

	public static double rootMeanSquare(double[] nums) {
	/*Calculates the root mean square of the input array.
	*
	* INPUT: Array to */
		double sum = 0.0;
		for (double num : nums)
			sum += num * num;
		return Math.sqrt(sum / nums.length);
	}

	public double get_std(double[] hrv){
		DescriptiveStatistics stats = new DescriptiveStatistics();
		for(int i=0;i<hrv.length;i++){
			stats.addValue(hrv[i]);
		}
		double res = stats.getStandardDeviation();
		return res;
	}


	public double[] get_diff(double[] signal){
	double[] diff_signal = new double[signal.length];

	for (int i=0;i<diff_signal.length-1;i++){
		diff_signal[i] = signal[i]-signal[i+1];
	}
		return diff_signal;
	}

	public double get_nn50(double[] signal){
		double nn50 = 0;
		int msInt = 50/1000;
		double[] over_msInt = new double[signal.length];

		for(int i=0;i<signal.length;i++){
			if(over_msInt[i] != 0 && over_msInt[i]>msInt){
				nn50 += over_msInt[i];
			}
		}
		return nn50;
	}



}
	// EMPATICA E4 stuff //TODO: Discuss wether to keep it in our not since we're not using it
	/*@Override
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
	}*/

	/*@Override
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


}*/

