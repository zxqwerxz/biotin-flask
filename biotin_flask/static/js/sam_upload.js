/**
 * @fileOverview Custom JQuery code loaded with html page.
 * @see 'sam/upload.html'
 * @author Jeffrey Zhou
 * @version 0.0.1
 */

 /*****************************************************************************
  * Main Method
  ****************************************************************************/

$(function() {

  var NO_FILES_UPLOADED = '<tr><td>No files currently uploaded.</td></tr>'

  // Delete a single SAM file when on badge click
  $('#samfileTable tbody span.link-hover').click(function(event) {
    var filename = $(this).attr('data-filename');
    data = { 'file' : filename };
    var that = $(this)
    $.ajax({
      type: 'POST',
      url: $SCRIPT_ROOT + '/sam/delete',
      data: JSON.stringify(data),
      contentType: 'application/json',
      success: function(result) {
        that.closest('tr').remove();
        if ($('#samfileTable tbody tr').length == 0) {
          $('#samfileTable tbody').append(NO_FILES_UPLOADED);
          $('#deleteAllButton').remove();
        }
      }
    });
  });

  // Delete all SAM files on button click
  $('#deleteAllButton').click(function(event) {
    $.ajax({
      type: 'POST',
      url: $SCRIPT_ROOT + '/sam/delete_allbam',
      data: JSON.stringify({ }),
      contentType: 'application/json',
      success: function(result) {
        $('#samfileTable tbody tr').remove();
        if ($('#samfileTable tbody tr').length == 0) {
          $('#samfileTable tbody').append(NO_FILES_UPLOADED);
          $('#deleteAllButton').remove();
        }
      }
    });
  });

  // Delete a single FASTA file when on badge click
  $('#fastafileTable tbody span.link-hover').click(function(event) {
    var filename = $(this).attr('data-filename');
    data = { 'file' : filename };
    var that = $(this)
    $.ajax({
      type: 'POST',
      url: $SCRIPT_ROOT + '/sam/delete',
      data: JSON.stringify(data),
      contentType: 'application/json',
      success: function(result) {
        that.closest('tr').remove();
        if ($('#fastafileTable tbody tr').length == 0) {
          $('#fastafileTable tbody').append(NO_FILES_UPLOADED);
          $('#fastaInput').parent().show();
        }
      }
    });
  });

  // Show the spinner where the submit button is clicked
  $('#submitButton').click(function(event) {
    $('#spinner').show();
    $('#page-cover').fadeIn(300);
  });
});
