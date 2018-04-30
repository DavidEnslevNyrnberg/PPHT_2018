package com.example.erik.acutestressapp;

import android.Manifest;
import android.app.Activity;
import android.bluetooth.BluetoothAdapter;
import android.bluetooth.BluetoothDevice;
import android.content.DialogInterface;
import android.content.Intent;
import android.content.pm.PackageManager;
import umich.cse.yctung.androidlibsvm.LibSVM;
import android.content.res.AssetFileDescriptor;
import android.graphics.Color;
import android.net.Uri;

import android.os.Bundle;


import android.support.v7.app.AppCompatActivity;

import android.view.View;
import android.widget.TextView;



import com.empatica.empalink.ConnectionNotAllowedException;
import com.empatica.empalink.EmpaDeviceManager;
import com.empatica.empalink.config.EmpaSensorStatus;
import com.empatica.empalink.config.EmpaSensorType;
import com.empatica.empalink.config.EmpaStatus;
import com.empatica.empalink.delegate.EmpaDataDelegate;
import com.empatica.empalink.delegate.EmpaStatusDelegate;
import com.jjoe64.graphview.GraphView;
import com.jjoe64.graphview.helper.DateAsXAxisLabelFormatter;
import com.jjoe64.graphview.series.DataPoint;
import com.jjoe64.graphview.series.LineGraphSeries;
import com.jjoe64.graphview.Viewport;
import com.opencsv.CSVReader;
import com.opencsv.bean.CsvToBeanBuilder;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;


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


	int lastX = 5;
	double eda_fs;
	double bvp_fs;
	private double[] edaData;
	private double[] bvpData;
	java.util.Date time;

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
		graph.getViewport().setScrollable(true); // enables horizontal scrolling
		graph.getViewport().setScrollableY(true); // enables vertical scrolling
		graph.getViewport().setScalable(true); // enables horizontal zooming and scrolling
		graph.getViewport().setScalableY(true); // enables vertical zooming and scrolling



		// Initialization of the grap series
		series_EDA = new LineGraphSeries<>();
		graph.addSeries(series_EDA);
		// styling series
		series_EDA.setColor(Color.argb(255,0,255,255));
		series_EDA.setDrawDataPoints(true);
		series_EDA.setThickness(3);
		/*graph.getGridLabelRenderer().setLabelFormatter(new DateAsXAxisLabelFormatter(getActivity()));
		graph.getGridLabelRenderer().setNumHorizontalLabels(3); // only 4 because of the space
		graph.getViewport().setMinX(getTime());*/
		// as we use dates as labels, the human rounding to nice readable numbers
		// is not necessary
		graph.getGridLabelRenderer().setHumanRounding(false);

		series_BVP = new LineGraphSeries<>();
		graph.addSeries(series_BVP);
		// styling series
		series_BVP.setColor(Color.argb(255,255,50,0));
		series_BVP.setDrawDataPoints(true);
		series_BVP.setThickness(3);

		edaData = getEDA();
		bvpData = getBVP();

	}



	@Override
	protected void onResume() {
		super.onResume();
		// Simulating real time with thread that append data to the graph
		new Thread(new Runnable() {

			@Override
			public void run() {
				// we add 10000 new entries
				for (int i = 0; i < 1000; i++) {
					runOnUiThread(new Runnable() {

						@Override
						public void run() {
							addEntry();
						}
					});


					// sleep to slow down the add of entries
					try {
						Thread.sleep(600);
					} catch (InterruptedException e) {
						// manage error ...
					}
				}
			}
		}).start();

	}

	// add data to graph
	private void addEntry() {
		// TODO: splitting EDA and BVP to add at correct sample rate
		// here, we choose to display max 100 points on the viewport and we scroll to end
		series_EDA.appendData(new DataPoint(lastX++, edaData[lastX++]), true, 100);
		series_BVP.appendData(new DataPoint(lastX++, bvpData[lastX++]), true, 100);

	}

	// Function to read the EDA file and do the pre-processing
	// TODO: Finish pre-processing
	private double[] getEDA(){
		// return array
		double[] res = new double[10000];

		int i = 0;

		try {

			// CSV read
			CSVReader reader = new CSVReader(new InputStreamReader(getAssets().open("eda.csv")));
			String[] next;
			while ((next = reader.readNext()) != null){
				res[i] = Double.parseDouble(next[0]);
				i = i + 1;
			}
			eda_fs = res[1];

		} catch (IOException e) {
			e.printStackTrace();
		}
		// LP filtering
		res = fourierLowPassFilter(res, 150, bvp_fs);
		return res;
	}
	// Function to read the BVP file and do the pre-processing
	// TODO: Finish pre-processing
	private double[] getBVP(){

		double[] res2 = new double[70000];
		int j = 0;

		try {

			// CSV read
			CSVReader reader = new CSVReader(new InputStreamReader(getAssets().open("bvp.csv")));
			String[] next;
			while ((next = reader.readNext()) != null){
				res2[j] = Double.parseDouble(next[0]);
				j = j + 1;
			}
			bvp_fs = res2[1];

		} catch (IOException e) {
			e.printStackTrace();
		}
		// LP filtering
		res2 = fourierLowPassFilter(res2, 150, bvp_fs);
		return res2;
	}
	// TODO: time converting
	private java.util.Date getTime(String filename){

		float[] time_x = new float[70000];


		try {

			// Reads the file and
			CSVReader reader = new CSVReader(new InputStreamReader(getAssets().open(filename)));
			String[] next;
			next = reader.readNext();
			time_x[0] = Float.parseFloat(next[0]);
			time = new java.util.Date((long)time_x[0]*1000);

		} catch (IOException e) {
			e.printStackTrace();
		}
		return time;
	}


	// Update text in Stress status
	public void updateTextView(String toThis) {
		TextView textView = statusLabel;
		textView.setText(toThis);
	}

	public void alarmHistory(View v){
		Intent i = new Intent(getApplicationContext(),alarmHistory.class);
		//StartActivity with the intent i
		startActivity(i);

	}



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



	/*public List<String[]> ppg2peak(){ // TODO: port ppg2peak from matlab
	}

	public List<String[]> getHRV(){ // TODO: port get_hrv from matlab

	}*/


// Moving average filter
public float[] movingAverageFilter(float[] list, int Arrsize){
	final int SMOOTH_FACTOR_MAA = 2;//increase for better results   but hits cpu bad

	int listSize = list.length; //input list
	int iterations = listSize / SMOOTH_FACTOR_MAA;
	float[] gList = new float[Arrsize];
	for (int i = 0, node = 0; i < iterations; i++) {
		float num = 0;
		for (int k = node; k < node + SMOOTH_FACTOR_MAA; k++) {
			num = num + list[k];
		}
		node = node + SMOOTH_FACTOR_MAA;
		num = num / SMOOTH_FACTOR_MAA;
		gList[i] = (num);//out put list
	}
	return gList;

}






}
































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

