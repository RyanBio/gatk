package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Train a Convolutional Neural Network (CNN) for filtering variants.
 * This tool expects requires training data generated by {@link CNNVariantWriteTensors}.
 *
 *
 * <h3>Inputs</h3>
 * <ul>
 *      <li>data-dir The training data created by {@link CNNVariantWriteTensors}.</li>
 *      <li>The tensor-name argument determines what types of tensors the model will expect.
 *      Set it to "reference" for 1D tensors or "read_tensor" for 2D tensors.</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 * <li>output-dir The model weights file and semantic configuration json are saved here.
 *  This default to the current working directory.</li>
 * <li>model-name The name for your model.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Train a 1D CNN on Reference Tensors</h4>
 * <pre>
 * gatk CNNVariantTrain \
 *   -tensor-type reference \
 *   -input-tensors-dir my_tensor_folder \
 *   -model-name my_1d_model
 * </pre>
 *
 * <h4>Train a 2D CNN on Read Tensors</h4>
 * <pre>
 * gatk CNNVariantTrain \
 *   -input-tensors-dir my_tensor_folder \
 *   -tensor-type read-tensor \
 *   -model-name my_2d_model
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Train a CNN model for filtering variants",
        oneLineSummary = "Train a CNN model for filtering variants",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class CNNVariantTrain extends CommandLineProgram {

    @Argument(fullName = "input-tensor-dir", shortName = "input-tensor-dir", doc = "Directory of training tensors to create.")
    private String inputTensorDir;

    @Argument(fullName = "output-dir", shortName = "output-dir", doc = "Directory where models will be saved, defaults to current working directory.", optional = true)
    private String outputDir = "./";

    @Argument(fullName = "tensor-type", shortName = "tensor-type", doc = "Name of the tensors to generate, reference for 1D reference tensors and read_tensor for 2D tensors.", optional = true)
    private TensorType tensorType = TensorType.reference;

    @Argument(fullName = "model-name", shortName = "model-name", doc = "Name of the model to be trained.", optional = true)
    private String modelName = "variant_filter_model";

    @Argument(fullName = "epochs", shortName = "epochs", doc = "Maximum number of training epochs.", optional = true, minValue = 0)
    private int epochs = 10;

    @Argument(fullName = "training-steps", shortName = "training-steps", doc = "Number of training steps per epoch.", optional = true, minValue = 0)
    private int trainingSteps = 10;

    @Argument(fullName = "validation-steps", shortName = "validation-steps", doc = "Number of validation steps per epoch.", optional = true, minValue = 0)
    private int validationSteps = 2;

    @Argument(fullName = "image-dir", shortName = "image-dir", doc = "Path where plots and figures are saved.", optional = true)
    private String imageDir;

    @Advanced
    @Argument(fullName = "channels-last", shortName = "channels-last", doc = "Store the channels in the last axis of tensors, tensorflow->true, theano->false", optional = true)
    private boolean channelsLast = true;

    @Advanced
    @Argument(fullName = "annotation-set", shortName = "annotation-set", doc = "Which set of annotations to use.", optional = true)
    private String annotationSet = "best_practices";

    // Start the Python executor. This does not actually start the Python process, but fails if python can't be located
    final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);


    @Override
    protected void onStartup() {
        PythonScriptExecutor.checkPythonEnvironmentForPackage("vqsr_cnn");
    }

    @Override
    protected Object doWork() {
        final Resource pythonScriptResource = new Resource("training.py", FilterVariantTranches.class);
        List<String> arguments = new ArrayList<>(Arrays.asList(
                "--data_dir", inputTensorDir,
                "--output_dir", outputDir,
                "--tensor_name", tensorType.name(),
                "--annotation_set", annotationSet,
                "--epochs", Integer.toString(epochs),
                "--training_steps", Integer.toString(trainingSteps),
                "--validation_steps", Integer.toString(validationSteps),
                "--id", modelName));

        if(channelsLast){
            arguments.add("--channels_last");
        } else {
            arguments.add("--channels_first");
        }

        if(imageDir != null){
            arguments.addAll(Arrays.asList("--image_dir", imageDir));
        }

        if (tensorType == TensorType.reference) {
            arguments.addAll(Arrays.asList("--mode", "train_on_reference_tensors_and_annotations"));
        } else if (tensorType == TensorType.read_tensor) {
            arguments.addAll(Arrays.asList("--mode", "train_small_model_on_read_tensors_and_annotations"));
        } else {
            throw new GATKException("Unknown tensor mapping mode:"+ tensorType.name());
        }

        logger.info("Args are:"+ Arrays.toString(arguments.toArray()));
        final boolean pythonReturnCode = pythonExecutor.executeScript(
                pythonScriptResource,
                null,
                arguments
        );
        return pythonReturnCode;
    }

}
